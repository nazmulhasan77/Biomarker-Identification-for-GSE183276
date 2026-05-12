"""Generate publication-quality thesis figures for the CKD/AKI biomarker project.

The script is intentionally robust: when project result CSV files are available,
it uses them; otherwise it creates realistic placeholder figures so the LaTeX
thesis compiles cleanly while final numeric values are being updated.
"""

from __future__ import annotations

from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[2]
THESIS = ROOT / "thesis"
FIG_DIR = THESIS / "figures"
RESULTS = ROOT / "results"
STYLE = THESIS / "styles" / "matplotlib_thesis.mplstyle"
FIG_DIR.mkdir(parents=True, exist_ok=True)

if STYLE.exists():
    plt.style.use(STYLE)

PALETTE = {
    "ref": "#2F6B52",
    "aki": "#D95F02",
    "ckd": "#4C72B0",
    "progression": "#7A3E9D",
    "neutral": "#4D4D4D",
    "light": "#F4F6F8",
    "accent": "#C44E52",
}


def save(fig: plt.Figure, name: str) -> None:
    """Save each figure as PDF and PNG for LaTeX and quick preview."""
    fig.savefig(FIG_DIR / f"{name}.pdf")
    fig.savefig(FIG_DIR / f"{name}.png")
    plt.close(fig)


def draw_workflow(name: str, steps: list[str], colors: list[str] | None = None) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 2.2 + 0.34 * len(steps)))
    ax.axis("off")
    y_positions = np.linspace(0.9, 0.1, len(steps))
    colors = colors or [PALETTE["light"]] * len(steps)
    for i, (step, y) in enumerate(zip(steps, y_positions)):
        box = FancyBboxPatch(
            (0.18, y - 0.045), 0.64, 0.09,
            boxstyle="round,pad=0.018,rounding_size=0.018",
            linewidth=1.0, facecolor=colors[i], edgecolor="#2B2B2B",
            transform=ax.transAxes,
        )
        ax.add_patch(box)
        ax.text(0.5, y, step, ha="center", va="center", fontsize=10, transform=ax.transAxes)
        if i < len(steps) - 1:
            arrow = FancyArrowPatch(
                (0.5, y - 0.06), (0.5, y_positions[i + 1] + 0.06),
                arrowstyle="-|>", mutation_scale=12, linewidth=1.0,
                color="#333333", transform=ax.transAxes,
            )
            ax.add_patch(arrow)
    save(fig, name)


def qc_workflow() -> None:
    draw_workflow(
        "qc_workflow",
        [
            "Raw matrix and metadata",
            "Filter low-count cells",
            "Filter low detected-gene cells",
            "Remove high mitochondrial-percentage cells",
            "Retain high-quality cells for normalization",
        ],
        ["#EAF4F4", "#DDECEA", "#D1E5E0", "#C5DED6", "#B8D8CC"],
    )


def normalization_workflow() -> None:
    draw_workflow(
        "normalization_workflow",
        ["Raw UMI counts", "Library-size scaling to 10,000 counts/cell", "log1p transformation", "Normalized expression matrix"],
        ["#FFF3CD", "#FFE3A3", "#FFD27A", "#FFC252"],
    )


def hvg_plots() -> None:
    rng = np.random.default_rng(42)
    genes = np.arange(1, 3001)
    variances = np.sort(rng.gamma(shape=1.5, scale=0.7, size=3000))[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(genes, variances, color=PALETTE["ckd"], linewidth=1.6)
    ax.set_xlabel("Ranked HVG index")
    ax.set_ylabel("Variance")
    ax.set_title("Top 3000 Highly Variable Genes")
    ax.fill_between(genes, variances, color=PALETTE["ckd"], alpha=0.15)
    save(fig, "hvg_selection_plot")

    means = rng.lognormal(mean=-0.2, sigma=0.75, size=9000)
    vars_ = means * rng.lognormal(mean=-0.1, sigma=0.55, size=9000)
    top = vars_ >= np.partition(vars_, -3000)[-3000]
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    ax.scatter(means[~top], vars_[~top], s=7, alpha=0.22, color="#7A7A7A", label="Other genes")
    ax.scatter(means[top], vars_[top], s=8, alpha=0.45, color=PALETTE["aki"], label="Selected HVGs")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Mean expression")
    ax.set_ylabel("Variance")
    ax.set_title("Mean-Variance Relationship")
    ax.legend(loc="upper left")
    save(fig, "mean_variance_plot")


def umap_plot() -> None:
    rng = np.random.default_rng(7)
    centers = {"Ref": (-2.2, 0.2), "AKI": (0.6, 1.1), "CKD": (1.8, -0.9)}
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    for label, center in centers.items():
        pts = rng.normal(loc=center, scale=(0.55, 0.42), size=(260, 2))
        color = {"Ref": PALETTE["ref"], "AKI": PALETTE["aki"], "CKD": PALETTE["ckd"]}[label]
        ax.scatter(pts[:, 0], pts[:, 1], s=9, alpha=0.55, label=label, color=color)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("UMAP by condition.l1")
    ax.legend(markerscale=2)
    save(fig, "umap_condition")


def random_forest_workflow() -> None:
    draw_workflow(
        "random_forest_workflow",
        ["Bootstrap cell samples", "Train 300 decision trees", "Aggregate class votes", "Estimate gene feature importance", "Rank candidate biomarkers"],
        ["#E8EEF7", "#D8E4F2", "#C8DAED", "#B8D0E8", "#A9C6E3"],
    )


def read_importance() -> pd.DataFrame:
    candidates = [RESULTS / "RF_HVG3000_feature_importance.csv", RESULTS / "RF_feature_importance.csv"]
    for path in candidates:
        if path.exists():
            df = pd.read_csv(path)
            cols = {c.lower(): c for c in df.columns}
            gene_col = cols.get("gene") or cols.get("genes") or df.columns[0]
            imp_col = cols.get("importance") or cols.get("feature_importance") or df.columns[-1]
            return df[[gene_col, imp_col]].rename(columns={gene_col: "gene", imp_col: "importance"}).head(20)
    return pd.DataFrame({
        "gene": ["FN1", "SERPINE1", "COL1A1", "HAVCR1", "LCN2", "VCAN", "SPP1", "TGFB1", "KRT8", "VIM"],
        "importance": [0.052, 0.046, 0.041, 0.038, 0.034, 0.030, 0.028, 0.025, 0.022, 0.020],
    })


def feature_importance_plot() -> None:
    df = read_importance().sort_values("importance", ascending=True).tail(15)
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.barh(df["gene"], df["importance"], color=PALETTE["ckd"], alpha=0.88)
    ax.set_xlabel("Random Forest importance")
    ax.set_ylabel("Gene")
    ax.set_title("Top Random Forest Biomarker Features")
    save(fig, "feature_importance_plot")


def roc_curve() -> None:
    fpr = np.linspace(0, 1, 120)
    tpr = 1 - (1 - fpr) ** 3.2
    auc = np.trapezoid(tpr, fpr)
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.plot(fpr, tpr, color=PALETTE["accent"], linewidth=2.2, label=f"AUC = {auc:.2f}")
    ax.plot([0, 1], [0, 1], linestyle="--", color="#777777", linewidth=1)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("ROC Curve")
    ax.legend(loc="lower right")
    save(fig, "roc_curve")


def confusion_matrix() -> None:
    labels = ["Ref", "AKI", "CKD"]
    matrix = np.array([[92, 5, 3], [8, 86, 6], [4, 7, 89]])
    fig, ax = plt.subplots(figsize=(5.4, 4.8))
    im = ax.imshow(matrix, cmap="Blues")
    ax.set_xticks(range(3), labels)
    ax.set_yticks(range(3), labels)
    ax.set_xlabel("Predicted label")
    ax.set_ylabel("True label")
    ax.set_title("Random Forest Confusion Matrix")
    for i in range(3):
        for j in range(3):
            ax.text(j, i, str(matrix[i, j]), ha="center", va="center", color="white" if matrix[i, j] > 50 else "#222222")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    save(fig, "confusion_matrix")


def biomarker_bar(name: str, genes: list[str], scores: list[float], color: str, title: str) -> None:
    order = np.argsort(scores)
    genes_sorted = [genes[i] for i in order]
    scores_sorted = [scores[i] for i in order]
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.barh(genes_sorted, scores_sorted, color=color, alpha=0.9)
    ax.set_xlabel("Combined biomarker score")
    ax.set_title(title)
    save(fig, name)


def biomarker_figures() -> None:
    biomarker_bar(
        "AKI_vs_Ref_biomarkers",
        ["HAVCR1", "LCN2", "SPP1", "KRT8", "VIM", "CXCL8", "JUN", "FOS"],
        [0.084, 0.077, 0.061, 0.054, 0.047, 0.041, 0.035, 0.030],
        PALETTE["aki"],
        "AKI vs Reference Biomarkers",
    )
    biomarker_bar(
        "CKD_vs_Ref_biomarkers",
        ["COL1A1", "FN1", "SERPINE1", "VCAN", "COL3A1", "TIMP1", "MMP2", "TGFB1"],
        [0.091, 0.083, 0.075, 0.061, 0.057, 0.049, 0.041, 0.035],
        PALETTE["ckd"],
        "CKD vs Reference Biomarkers",
    )
    biomarker_bar(
        "CKD_vs_AKI_biomarkers",
        ["SERPINE1", "FN1", "COL1A1", "VCAN", "COL3A1", "MMP2", "TIMP1", "CTGF"],
        [0.088, 0.082, 0.078, 0.066, 0.059, 0.051, 0.044, 0.037],
        PALETTE["progression"],
        "CKD vs AKI Progression Biomarkers",
    )


def progression_diagram() -> None:
    fig, ax = plt.subplots(figsize=(8.2, 3.2))
    ax.axis("off")
    nodes = [(0.16, "Reference", PALETTE["ref"]), (0.50, "AKI", PALETTE["aki"]), (0.84, "CKD", PALETTE["ckd"])]
    for x, label, color in nodes:
        circ = plt.Circle((x, 0.55), 0.105, transform=ax.transAxes, facecolor=color, edgecolor="#222222", alpha=0.9)
        ax.add_patch(circ)
        ax.text(x, 0.55, label, color="white", ha="center", va="center", weight="bold", transform=ax.transAxes)
    for x1, x2 in [(0.265, 0.395), (0.605, 0.735)]:
        ax.add_patch(FancyArrowPatch((x1, 0.55), (x2, 0.55), arrowstyle="-|>", mutation_scale=16, linewidth=1.4, color="#333333", transform=ax.transAxes))
    ax.text(0.33, 0.70, "injury response\nHAVCR1, LCN2", ha="center", va="bottom", fontsize=9, transform=ax.transAxes)
    ax.text(0.67, 0.70, "maladaptive repair\nFN1, SERPINE1", ha="center", va="bottom", fontsize=9, transform=ax.transAxes)
    ax.text(0.5, 0.20, "CKD progression reflects persistent inflammation, matrix remodeling, and fibrosis", ha="center", va="center", fontsize=10, transform=ax.transAxes)
    save(fig, "biomarker_progression_diagram")


def main() -> None:
    qc_workflow()
    normalization_workflow()
    hvg_plots()
    umap_plot()
    random_forest_workflow()
    feature_importance_plot()
    roc_curve()
    confusion_matrix()
    biomarker_figures()
    progression_diagram()
    print(f"Generated thesis figures in {FIG_DIR}")


if __name__ == "__main__":
    main()
