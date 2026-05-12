"""GSE183276 single-cell RNA-seq CKD biomarker pipeline.

Run from the project root:
    python biomarker_pipeline.py
"""

from pathlib import Path
import sys
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder


DATA_DIR = Path("data")
RESULTS_DIR = Path("results")

COUNTS_FILE = DATA_DIR / "counts.mtx"
GENES_FILE = DATA_DIR / "genes.tsv"
BARCODES_FILE = DATA_DIR / "barcodes.tsv"
METADATA_FILE = (
    DATA_DIR / "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_03282022.txt"
)

LABEL_COLUMN = "condition.l1"
CELLTYPE_COLUMN = "subclass.l1"
OPTIONAL_COLUMNS = ["class", "state", "structure", "condition.l2", "condition.l3"]
UMAP_COLUMNS = [LABEL_COLUMN, CELLTYPE_COLUMN, "class", "state"]


def require_files(paths):
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required input file(s):\n"
            + "\n".join(f"  - {path}" for path in missing)
            + "\nRun this script from the project root and keep all inputs in data/."
        )


def read_single_column_tsv(path, name):
    values = pd.read_csv(path, sep="\t", header=None, dtype=str).iloc[:, 0]
    values = values.fillna("").astype(str)
    if values.empty:
        raise ValueError(f"{path} is empty; cannot load {name}.")
    return values


def make_unique(values):
    """Return unique string labels while preserving the first occurrence unchanged."""
    counts = {}
    unique = []
    for value in values:
        value = str(value)
        if value not in counts:
            counts[value] = 0
            unique.append(value)
        else:
            counts[value] += 1
            unique.append(f"{value}-{counts[value]}")
    return unique


def load_counts():
    print("Loading raw count matrix...")
    genes = read_single_column_tsv(GENES_FILE, "genes")
    barcodes = read_single_column_tsv(BARCODES_FILE, "barcodes")

    counts = mmread(COUNTS_FILE)
    if not sparse.issparse(counts):
        counts = sparse.coo_matrix(counts)

    # Matrix Market files from 10X-style data are commonly genes x cells.
    # Transpose to cells x genes as requested.
    counts = counts.tocsr().T.tocsr()

    expected_shape = (len(barcodes), len(genes))
    if counts.shape != expected_shape:
        raise ValueError(
            "Count matrix dimensions do not match barcodes/genes after transpose.\n"
            f"  counts.mtx after transpose: {counts.shape}\n"
            f"  expected cells x genes:     {expected_shape}\n"
            "Check that counts.mtx is genes x cells and that genes.tsv/barcodes.tsv match it."
        )

    adata = ad.AnnData(X=counts)
    adata.obs_names = make_unique(barcodes)
    adata.var_names = make_unique(genes)
    adata.var["gene_name"] = adata.var_names

    print(f"Loaded AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    return adata


def load_metadata(barcodes):
    print("Loading metadata...")
    metadata = pd.read_csv(METADATA_FILE, sep="\t", index_col=0, low_memory=False)
    metadata.index = metadata.index.astype(str)

    common = pd.Index(barcodes).intersection(metadata.index)
    if len(common) == 0:
        raw_metadata = pd.read_csv(METADATA_FILE, sep="\t", low_memory=False)
        barcode_set = set(map(str, barcodes))

        best_column = None
        best_overlap = 0
        for column in raw_metadata.columns:
            overlap = raw_metadata[column].astype(str).isin(barcode_set).sum()
            if overlap > best_overlap:
                best_column = column
                best_overlap = overlap

        if best_column is None or best_overlap == 0:
            raise ValueError(
                "Could not match metadata rows to cell barcodes. "
                "Expected the first metadata column, or another column, to contain barcodes."
            )

        metadata = raw_metadata.set_index(best_column)
        metadata.index = metadata.index.astype(str)
        common = pd.Index(barcodes).intersection(metadata.index)

    print(f"Metadata rows: {metadata.shape[0]:,}; matched cells: {len(common):,}")
    return metadata


def attach_metadata(adata, metadata):
    common_cells = adata.obs_names.intersection(metadata.index)
    if len(common_cells) == 0:
        raise ValueError("No common cell barcodes found between AnnData and metadata.")

    adata = adata[common_cells].copy()
    metadata = metadata.loc[common_cells].copy()

    columns_to_keep = [
        column
        for column in [LABEL_COLUMN, CELLTYPE_COLUMN, *OPTIONAL_COLUMNS]
        if column in metadata.columns
    ]
    other_columns = [column for column in metadata.columns if column not in columns_to_keep]
    metadata = metadata[columns_to_keep + other_columns]

    adata.obs = adata.obs.join(metadata, how="left")

    for required_column in [LABEL_COLUMN, CELLTYPE_COLUMN]:
        if required_column not in adata.obs.columns:
            raise KeyError(
                f"Required metadata column '{required_column}' was not found in {METADATA_FILE}."
            )

    print(f"After metadata join: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print_label_counts(adata, LABEL_COLUMN)
    print_label_counts(adata, CELLTYPE_COLUMN)
    return adata


def print_label_counts(adata, column):
    if column in adata.obs.columns:
        print(f"\n{column} counts:")
        print(adata.obs[column].value_counts(dropna=False).to_string())


def run_quality_control(adata):
    print("\nRunning quality control...")
    if "percent.mt" in adata.obs.columns:
        adata.obs["percent.mt"] = pd.to_numeric(adata.obs["percent.mt"], errors="coerce")
        print("Using existing metadata column: percent.mt")
    else:
        print("percent.mt not found; calculating from genes starting with 'MT-'.")
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        if adata.var["mt"].sum() == 0:
            warnings.warn(
                "No mitochondrial genes starting with 'MT-' were found; percent.mt will be 0."
            )
            adata.obs["percent.mt"] = 0.0
        else:
            sc.pp.calculate_qc_metrics(
                adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
            )
            adata.obs["percent.mt"] = adata.obs["pct_counts_mt"]

    before = adata.shape
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["percent.mt"].fillna(np.inf) < 20].copy()

    print(f"QC before: {before[0]:,} cells x {before[1]:,} genes")
    print(f"QC after:  {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print_label_counts(adata, LABEL_COLUMN)
    return adata


def normalize_and_select_hvgs(adata):
    print("\nNormalizing and selecting highly variable genes...")
    adata.layers["raw_counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["log_normalized"] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    hvg_count = int(adata.var["highly_variable"].sum())
    if hvg_count == 0:
        raise ValueError("No highly variable genes were selected.")

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    print(f"Selected HVGs: {adata_hvg.n_vars:,}")
    return adata, adata_hvg


def run_visualization(adata_hvg):
    print("\nRunning PCA, neighbors, and UMAP...")
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, svd_solver="arpack")
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata_hvg)

    for column in UMAP_COLUMNS:
        if column in adata_hvg.obs.columns:
            output_file = RESULTS_DIR / f"UMAP_{column.replace('.', '_')}.png"
            sc.pl.umap(
                adata_hvg,
                color=column,
                show=False,
                save=None,
                frameon=False,
            )
            import matplotlib.pyplot as plt

            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            plt.close()
            print(f"Saved {output_file}")
        else:
            print(f"Skipping UMAP for missing column: {column}")

    return adata_hvg


def run_deg_analysis(adata):
    print("\nRunning differential expression analysis...")
    valid = adata.obs[LABEL_COLUMN].notna()
    adata_deg = adata[valid].copy()
    adata_deg.obs[LABEL_COLUMN] = adata_deg.obs[LABEL_COLUMN].astype(str).astype("category")

    if adata_deg.obs[LABEL_COLUMN].nunique() < 2:
        raise ValueError(
            f"DEG analysis requires at least two non-missing groups in '{LABEL_COLUMN}'."
        )

    sc.tl.rank_genes_groups(
        adata_deg,
        groupby=LABEL_COLUMN,
        method="wilcoxon",
        use_raw=False,
        pts=True,
    )
    deg_results = sc.get.rank_genes_groups_df(adata_deg, group=None)
    deg_results.to_csv(RESULTS_DIR / "DEG_results.csv", index=False)
    print(f"Saved DEG results: {deg_results.shape[0]:,} rows")
    return deg_results


def prepare_feature_matrix(adata_hvg):
    valid = adata_hvg.obs[LABEL_COLUMN].notna()
    ml_data = adata_hvg[valid].copy()

    y_labels = ml_data.obs[LABEL_COLUMN].astype(str).values
    if pd.Series(y_labels).nunique() < 2:
        raise ValueError(
            f"Random Forest requires at least two classes in '{LABEL_COLUMN}'."
        )

    X = ml_data.X
    if sparse.issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)

    return X, y_labels, ml_data.var_names.to_numpy()


def train_random_forest(adata_hvg):
    print("\nTraining Random Forest classifier...")
    X, y_labels, gene_names = prepare_feature_matrix(adata_hvg)

    label_encoder = LabelEncoder()
    y = label_encoder.fit_transform(y_labels)

    class_counts = pd.Series(y).value_counts()
    stratify = y if class_counts.min() >= 2 else None
    if stratify is None:
        warnings.warn(
            "At least one class has fewer than 2 cells; train/test split will not be stratified."
        )

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.2,
        random_state=42,
        stratify=stratify,
    )

    rf = RandomForestClassifier(
        n_estimators=300,
        random_state=42,
        n_jobs=-1,
        class_weight="balanced",
    )
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)

    report = classification_report(
        y_test,
        y_pred,
        target_names=label_encoder.classes_,
        zero_division=0,
    )
    (RESULTS_DIR / "classification_report.txt").write_text(report, encoding="utf-8")

    cm = pd.DataFrame(
        confusion_matrix(y_test, y_pred),
        index=[f"true_{label}" for label in label_encoder.classes_],
        columns=[f"pred_{label}" for label in label_encoder.classes_],
    )
    cm.to_csv(RESULTS_DIR / "confusion_matrix.csv")

    feature_importance = pd.DataFrame(
        {
            "gene": gene_names,
            "importance": rf.feature_importances_,
        }
    ).sort_values("importance", ascending=False)

    feature_importance.to_csv(
        RESULTS_DIR / "RandomForest_feature_importance.csv", index=False
    )
    feature_importance.head(100).to_csv(
        RESULTS_DIR / "Top_Biomarker_Genes_RF.csv", index=False
    )

    print("Saved classification report, confusion matrix, and RF feature importance.")
    return feature_importance


def build_final_biomarker_ranking(feature_importance, deg_results):
    print("\nBuilding final biomarker ranking...")
    deg = deg_results.copy()
    rename_map = {
        "names": "gene",
        "logfoldchanges": "logfoldchange",
        "pvals_adj": "adj_p_value",
    }
    deg = deg.rename(columns=rename_map)

    required_columns = {"gene", "logfoldchange", "adj_p_value"}
    missing = required_columns - set(deg.columns)
    if missing:
        raise KeyError(
            "DEG results are missing required column(s) for final ranking: "
            + ", ".join(sorted(missing))
        )

    deg["logfoldchange"] = pd.to_numeric(deg["logfoldchange"], errors="coerce")
    deg["adj_p_value"] = pd.to_numeric(deg["adj_p_value"], errors="coerce")

    merged = feature_importance.merge(deg, on="gene", how="inner")
    merged["combined_score"] = (
        merged["importance"]
        * merged["logfoldchange"].abs()
        * (-np.log10(merged["adj_p_value"].fillna(1.0) + 1e-300))
    )
    merged = merged.sort_values("combined_score", ascending=False)
    merged.to_csv(RESULTS_DIR / "Final_Biomarker_Candidates_RF_DEG.csv", index=False)
    print(f"Saved final biomarker ranking: {merged.shape[0]:,} rows")
    return merged


def save_processed_data(adata, adata_hvg):
    print("\nSaving processed AnnData objects...")
    adata.write_h5ad(RESULTS_DIR / "GSE183276_processed_all_genes.h5ad")
    adata_hvg.write_h5ad(RESULTS_DIR / "GSE183276_processed_HVG.h5ad")
    print("Saved processed .h5ad files.")


def main():
    sc.settings.verbosity = 2
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    require_files([COUNTS_FILE, GENES_FILE, BARCODES_FILE, METADATA_FILE])

    adata = load_counts()
    metadata = load_metadata(adata.obs_names)
    adata = attach_metadata(adata, metadata)
    adata = run_quality_control(adata)
    adata, adata_hvg = normalize_and_select_hvgs(adata)
    adata_hvg = run_visualization(adata_hvg)
    deg_results = run_deg_analysis(adata)
    feature_importance = train_random_forest(adata_hvg)
    build_final_biomarker_ranking(feature_importance, deg_results)
    save_processed_data(adata, adata_hvg)

    print("\nPipeline complete. Outputs are in the results/ folder.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        sys.exit(1)
