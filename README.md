# CKD Biomarker Pipeline for GSE183276

This project identifies candidate CKD-related biomarkers from the GSE183276 single-cell RNA-seq dataset using raw Matrix Market counts and metadata.

## How to run

```bash
pip install -r requirements.txt
python biomarker_pipeline.py
```

## Input files

Place these files in `data/`:

- `counts.mtx`
- `genes.tsv`
- `barcodes.tsv`
- `GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_03282022.txt`

The pipeline uses `condition.l1` as the disease/control label and `subclass.l1` as the cell type annotation. It also keeps `class`, `state`, `structure`, `condition.l2`, and `condition.l3` when available.

## Pipeline steps

1. Load raw counts, genes, and barcodes into an AnnData object.
2. Load and join metadata by cell barcode.
3. Run quality control with `min_genes >= 200`, `min_cells >= 3`, and `percent.mt < 20`.
4. Store raw counts, normalize to 10,000 counts per cell, and apply `log1p`.
5. Select the top 2,000 highly variable genes.
6. Run scaling, PCA, neighbors, and UMAP on HVGs.
7. Run Wilcoxon differential expression by `condition.l1`.
8. Train a balanced Random Forest classifier using HVG expression.
9. Merge Random Forest importance with DEG statistics for final biomarker ranking.
10. Save processed AnnData files.

## Output files

All outputs are written to `results/`:

- `DEG_results.csv`
- `classification_report.txt`
- `confusion_matrix.csv`
- `RandomForest_feature_importance.csv`
- `Top_Biomarker_Genes_RF.csv`
- `Final_Biomarker_Candidates_RF_DEG.csv`
- `UMAP_condition_l1.png`
- `UMAP_subclass_l1.png`
- `UMAP_class.png` if available
- `UMAP_state.png` if available
- `GSE183276_processed_all_genes.h5ad`
- `GSE183276_processed_HVG.h5ad`
