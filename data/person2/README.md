# Person 2. Preprocessing and quality control

The executed [Stage 3 notebook](../../notebooks/02_preprocessing.ipynb) reads the supplied GSE15852 expression and metadata CSVs, records their SHA-256 fingerprints and writes this evidence package. The [GEO study record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15852) describes the Affymetrix GPL96 assay and paired design.

## Requirements and evidence

| Requirement | Result | Evidence |
| --- | --- | --- |
| Align expression and metadata | All 86 sample IDs matched. Labels and title-derived patient IDs agreed. All 43 patient pairs were complete. | [Aligned metadata](clean_metadata.csv), executed notebook |
| Check missing, nonnumeric and duplicate values | Zero missing or nonnumeric expression cells, duplicate sample IDs, exact duplicate expression rows or duplicate sample profiles. | [Before/after checks](qc_before_after.csv) |
| Filter and transform | 21,225 input rows became 13,237 unique gene-symbol rows. The per-sample median combined repeated-symbol rows, reducing the row count by 7,988. The placeholder, exact-row and zero-variance filters removed zero rows. The export uses `log2(intensity + 1)`. | [Change counts](qc_changes.csv), [input-row audit](gene_row_audit.csv), [clean expression matrix](clean_expression.csv) |
| Screen unusual samples and document exclusions | Five samples crossed the lower correlation fence. None crossed the median intensity fences. All 86 remained. GSM398100 had the lowest median correlation (0.676). | [Sample QC table](sample_qc.csv), [QC plot](sample_qc.png) |
| Preserve raw counts where applicable | The input contains processed microarray signal intensities. The supplied input files remain available for comparison. | Notebook provenance and [dataset manifest](../dataset_manifest.md) |

The five review flags were GSM398074, GSM398100, GSM398123, GSM398154 and GSM398155. Both samples from Patient_0191 were flagged. The QC plot also shows two bands of median input intensity, each containing both tissue classes. The supplied matrix lacks raw chip diagnostics and batch labels, so the causes of these patterns remain unresolved.

## Handover to Persons 3, 4 and 6

- The expression CSV has gene symbols as rows and GEO sample IDs as columns. The metadata CSV has one row per sample in the same order, with `patient_id`, `condition` and `label` preserved. The symbols come from the acquisition export and include legacy labels.
- Person 3 can use the logged matrix directly for EDA. Person 4 can use it for exploratory paired analysis on the logged microarray scale. Match samples by `patient_id` and account for the QC flags in interpretation.
- Person 6 should group splits by patient. Fit any additional feature selection or scaling on training data only. The shared transformation is element-wise and uses a fixed formula. The cohort-wide constant-gene check removed zero features.
- Person 1 can integrate the Stage 3 notebook cells. For Colab, provide the two original CSV inputs to the runtime or reconstruct them through Stage 2. Each review output can be regenerated from those files.

