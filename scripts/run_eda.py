"""Generate the Group 10 Stage 4 figures from the Person 2 handover."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "person2"
OUTPUT = ROOT / "data" / "person3"
COLORS = {"Normal": "#146E79", "Tumor": "#AB462D"}
INK = "#202C3A"
GRID = "#E4E9EE"
plt.rcParams.update({"text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK, "axes.titleweight": "bold"})


def load_inputs():
    expression_path = INPUT / "clean_expression.csv"
    metadata_path = INPUT / "clean_metadata.csv"
    if not expression_path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError("Place clean_expression.csv and clean_metadata.csv in data/person2.")

    expression = pd.read_csv(expression_path, index_col="Gene_Symbol", keep_default_na=False)
    metadata = pd.read_csv(metadata_path).set_index("sample_id")
    if not expression.index.is_unique or not metadata.index.is_unique:
        raise ValueError("Expected unique gene symbols and sample IDs.")
    if set(expression.columns) != set(metadata.index):
        raise ValueError("Expression column IDs and metadata sample IDs must match.")
    if not np.isfinite(expression.to_numpy(dtype=float)).all():
        raise ValueError("Clean expression contains missing or invalid values.")
    metadata = metadata.reindex(expression.columns)
    if set(metadata["condition"]) != set(COLORS) or metadata["patient_id"].isna().any():
        raise ValueError("Expected Normal and Tumor labels and a patient ID for every sample.")
    pairs = pd.crosstab(metadata["patient_id"], metadata["condition"])
    if not pairs.eq(1).all().all():
        raise ValueError("Expected one Normal and one Tumor sample per patient.")
    return expression, metadata


def save_figure(figure, filename):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    path = OUTPUT / filename
    figure.savefig(path, dpi=160, facecolor="white")
    plt.close(figure)
    return path


def style_axis(ax, grid_axis="y"):
    ax.set_axisbelow(True)
    ax.grid(axis=grid_axis, color=GRID, linewidth=.8)
    ax.spines[["right", "top"]].set_visible(False)


def matched_samples(metadata, condition):
    return metadata.loc[metadata["condition"].eq(condition)].sort_values("patient_id").index


def distribution_plot(expression, metadata):
    fig, (ax, paired) = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained",
                                     gridspec_kw={"width_ratios": [1.4, 1]})
    bins = np.linspace(expression.min().min(), expression.max().max(), 65)
    for condition, color in COLORS.items():
        samples = matched_samples(metadata, condition)
        values = expression.loc[:, samples]
        ax.hist(values.to_numpy().ravel(), bins=bins, density=True, histtype="step",
                color=color, linewidth=2, linestyle="-" if condition == "Normal" else "--",
                label=f"{condition} (n={len(samples)})")
        medians = values.median(axis=0)
        print(f"{condition} per-sample median range: {medians.min():.2f}–{medians.max():.2f}")
    ax.set(title="A  Pooled expression distributions",
           xlabel="log₂(microarray signal intensity + 1)",
           ylabel="Density of gene–sample measurements")
    ax.legend(frameon=False)
    patient_medians = [expression.loc[:, matched_samples(metadata, tissue)].median().to_numpy()
                       for tissue in COLORS]
    paired.plot([0, 1], np.stack(patient_medians), color="#77838F", linewidth=.8, alpha=.46)
    for position, (tissue, color) in enumerate(COLORS.items()):
        paired.scatter(np.full(len(patient_medians[position]), position),
                       patient_medians[position], s=16, color=color, alpha=.8, zorder=3)
    paired.set(xlim=(-.24, 1.24), title="B  Sample medians by patient",
               ylabel="Median log₂(intensity + 1)",
               xlabel=f"Each line links one patient (n={len(patient_medians[0])})")
    paired.set_xticks([0, 1], list(COLORS))
    for panel in (ax, paired):
        style_axis(panel)
    fig.suptitle("GSE15852 expression ranges in matched breast tissue", fontsize=15, weight="bold")
    return save_figure(fig, "01_expression_distribution.png")


def pca_plot(expression, metadata):
    values = expression.to_numpy(dtype=float).T
    centered = values - values.mean(axis=0)
    # The sample Gram matrix yields the same PCA sample scores as SVD.
    eigenvalues, eigenvectors = np.linalg.eigh(centered @ centered.T)
    order = np.argsort(eigenvalues)[::-1][:2]
    variance = np.clip(eigenvalues, 0, None)
    explained = 100 * variance[order] / variance.sum()
    scores = eigenvectors[:, order] * np.sqrt(variance[order])
    for pc in range(scores.shape[1]):
        # Fix the arbitrary sign using the strongest gene loading.
        loadings = centered.T @ scores[:, pc]
        if loadings[np.argmax(np.abs(loadings))] < 0:
            scores[:, pc] *= -1

    fig, (ax, scree) = plt.subplots(1, 2, figsize=(11, 4.8), layout="constrained",
                                    gridspec_kw={"width_ratios": [2, 1]})
    for condition, color in COLORS.items():
        mask = metadata["condition"].eq(condition).to_numpy()
        ax.scatter(*scores[mask].T, color=color, alpha=.82, s=42,
                   label=f"{condition} (n={int(mask.sum())})")
    qc_file = INPUT / "sample_qc.csv"
    flagged = pd.Index([])
    if qc_file.is_file():
        qc = pd.read_csv(qc_file).set_index("sample_id")
        review = qc["review_flag"].astype(str).str.lower().eq("true")
        flagged = metadata.index.intersection(qc.index[review])
        mask = metadata.index.isin(flagged)
        if mask.any():
            ax.scatter(*scores[mask].T, facecolors="none", edgecolors="#27272A",
                       s=100, linewidths=1.2, label="Stage 3 QC flag")
    if len(flagged):
        flagged_positions = np.flatnonzero(metadata.index.isin(flagged))
        label_pos = flagged_positions[np.linalg.norm(scores[flagged_positions], axis=1).argmax()]
        ax.annotate(metadata.index[label_pos], scores[label_pos], xytext=(6, 7),
                    textcoords="offset points", fontsize=8)
    ax.set(title="A  Sample scores",
           xlabel=f"PC1 ({explained[0]:.1f}% variance)",
           ylabel=f"PC2 ({explained[1]:.1f}% variance)")
    ax.margins(.13)
    ax.legend(frameon=False)
    components = np.arange(1, min(10, len(variance)) + 1)
    percent = 100 * np.sort(variance)[::-1][:len(components)] / variance.sum()
    scree.bar(components, percent, color=["#365C74" if i <= 2 else "#BAC8D1" for i in components])
    scree.set(title="B  Variance by component", xlabel="Principal component",
              ylabel="Explained variance (%)", xticks=components)
    for panel in (ax, scree):
        style_axis(panel)
    fig.suptitle("PCA of GSE15852 expression profiles", fontsize=15, weight="bold")
    print(f"Variance explained: PC1 {explained[0]:.2f}%, PC2 {explained[1]:.2f}%.")
    positions = np.linalg.norm(scores, axis=1).argsort()[-3:][::-1]
    print("Furthest from the origin in PC1–PC2:", ", ".join(metadata.index[positions]))
    print(f"Stage 3 flags displayed: {len(flagged)}.")
    return save_figure(fig, "02_pca_samples.png")


def heatmap_plot(expression, metadata):
    variances = expression.var(axis=1)
    genes = variances.loc[variances.gt(0)].nlargest(30).index
    if genes.empty:
        raise ValueError("The expression matrix has no variable genes.")
    ordered = metadata.sort_values(["condition", "patient_id"])
    values = expression.loc[genes, ordered.index]
    scaled = values.sub(values.mean(axis=1), axis=0).div(values.std(axis=1), axis=0)
    normal_count = int(ordered["condition"].eq("Normal").sum())
    band = ordered["condition"].eq("Tumor").to_numpy(dtype=int)[None, :]

    fig = plt.figure(figsize=(12.2, 8.2), layout="constrained")
    axes = fig.add_gridspec(2, 3, height_ratios=[.4, 9], width_ratios=[22, 1, 6])
    annotation = fig.add_subplot(axes[0, 0])
    ax = fig.add_subplot(axes[1, 0])
    color_axis = fig.add_subplot(axes[1, 1])
    differences_axis = fig.add_subplot(axes[1, 2], sharey=ax)
    annotation.imshow(band, aspect="auto", cmap=ListedColormap(list(COLORS.values())),
                      vmin=0, vmax=1, interpolation="nearest")
    annotation.set_axis_off()
    annotation.legend(handles=[Patch(color=color, label=condition) for condition, color in COLORS.items()],
                      loc="lower center", bbox_to_anchor=(.5, .9), ncol=2, frameon=False)
    heat = ax.imshow(scaled, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                     interpolation="nearest")
    ax.axvline(normal_count - .5, color="#222222", linewidth=1)
    ax.set_yticks(range(len(genes)), genes, fontsize=8)
    ax.set_xticks([])
    ax.set(xlabel="Samples ordered by tissue, then patient ID",
           ylabel="Gene symbols ranked by variance")
    fig.suptitle(f"GSE15852 expression of the {len(genes)} most variable genes",
                 fontsize=15, weight="bold")
    fig.colorbar(heat, cax=color_axis, label="Within-gene z-score", extend="both")
    tumour = expression.loc[genes, matched_samples(metadata, "Tumor")].to_numpy()
    normal = expression.loc[genes, matched_samples(metadata, "Normal")].to_numpy()
    paired_change = np.median(tumour - normal, axis=1)
    differences_axis.barh(np.arange(len(genes)), paired_change,
                          color=np.where(paired_change >= 0, COLORS["Tumor"], COLORS["Normal"]),
                          height=.68)
    differences_axis.axvline(0, color=INK, linewidth=.8)
    differences_axis.tick_params(axis="y", left=False, labelleft=False)
    differences_axis.set(title="Matched pairs", xlabel="Median paired difference\nTumour − normal (log₂ scale)")
    style_axis(differences_axis, grid_axis="x")
    by_gene = pd.Series(paired_change, index=genes)
    print("Largest positive median paired differences:", ", ".join(by_gene.nlargest(3).index))
    print("Largest negative median paired differences:", ", ".join(by_gene.nsmallest(3).index))
    print("Five highest-variance genes:", ", ".join(genes[:5]))
    return save_figure(fig, "03_top_variable_genes_heatmap.png")


def main():
    expression, metadata = load_inputs()
    print(f"{expression.shape[0]:,} genes, {expression.shape[1]} samples, "
          f"{metadata['patient_id'].nunique()} matched pairs.")
    for plot in (distribution_plot, pca_plot, heatmap_plot):
        print(f"Saved {plot(expression, metadata).relative_to(ROOT)}")


if __name__ == "__main__":
    main()
