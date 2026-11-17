Accession & URL: NCBI GEO GSE15852 ([https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15852))

Download Date: November 17, 2026

Assay Type & Units: Microarray gene expression profiling (Affymetrix Human Genome U133A Array, GPL96); normalized signal intensity values (non-integer continuous expression values, not raw RNA-seq integer counts).   
Bioinformatics_Group_10_Six_Person_Work_Plan.pdf

Sample Counts & Pairing Status: 86 total samples consisting of 43 primary breast carcinoma samples (Tumor, label = 1) and 43 patient-matched adjacent normal breast tissue samples (Normal, label = 0) from 43 patients. Genuinely patient-matched (patient_id recorded in metadata.csv).
File Descriptions & Conventions:

expression_matrix.csv: Rows = Gene_Symbol (HGNC gene symbols), Columns = sample_id (GEO GSM accession IDs).   
Bioinformatics_Group_10_Six_Person_Work_Plan.pdf

metadata.csv: Columns = sample_id, sample_title, patient_id, condition (Tumor/Normal), label (1/0), source_tissue.

Suitability & Method Note for Teammates: Directly answers the tumor vs. adjacent normal breast tissue research question. Because values are normalized microarray intensities rather than raw RNA-seq counts, Person 4 must use a normalized/exploratory differential expression test (e.g., paired/two-sample t-test or Mann-Whitney U with Benjamini-Hochberg FDR) and explicitly note this assay limitation.
