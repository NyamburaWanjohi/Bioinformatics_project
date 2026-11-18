# Stage 4. Exploratory data analysis

The [executed notebook](../../notebooks/03_eda.ipynb) explores 13,237 logged gene expression features from 86 GSE15852 microarrays. There are 43 normal and 43 tumour samples from 43 matched patients. The notebook reads Person 2's cleaned expression and metadata together and uses the same plotting functions as the command below.

## Figures and interpretation

1. [Expression distributions](01_expression_distribution.png). The pooled distributions of log₂(intensity + 1) nearly overlap. A second panel links each patient's tumour and normal medians and shows how the medians vary within and between patients. The pooled view can conceal changes in individual genes.
2. [Principal component analysis](02_pca_samples.png). PC1 explains 33.85% of the variance and separates the samples into two broad bands that each contain both tissue classes. The variance panel shows how much the first ten components explain. GSM398074 is the labelled QC flag farthest from the origin. All five Stage 3 flags were retained. Batch information would be needed to assess whether the bands have a technical cause.
3. [Top-variable-gene heatmap](03_top_variable_genes_heatmap.png). The heatmap shows within-gene z-scores. The side panel shows the median of the 43 within-patient tumour-minus-normal differences for each gene on the logged scale. Among these 30 genes, KRT19 and AGR2 have large positive differences, while RBP4 and PPP1R1A have large negative differences. These are descriptive summaries.

## Run on Windows

1. Open the repository folder in VS Code and choose **Terminal → New Terminal**.
2. Check that `data/person2/` contains `clean_expression.csv` and `clean_metadata.csv`. Keeping `sample_qc.csv` in the same folder also reproduces the Stage 3 outlier markers.
3. With Python 3.12 installed, run the following from the repository root. If your Windows installation exposes Python as `python`, use `python` in place of `py -3.12`.

   ```powershell
   py -3.12 -m pip install -r requirements-eda.txt
   py -3.12 scripts/run_eda.py
   ```

4. Open the three PNGs in `data/person3/`. To rerun the notebook in VS Code, select a Python environment with these packages and a Jupyter kernel.

For the group notebook, place these three figures and their interpretations under Stage 4, or call the plot functions from the Stage 4 code. The heatmap's variance ranking is for exploration. Any feature selection used for modelling should be fitted on training data.
