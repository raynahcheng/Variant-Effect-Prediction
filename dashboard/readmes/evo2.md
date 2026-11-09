## About Evo2

Evo2 is a biological foundation model trained on DNA sequences from a wide range of organisms. Unlike traditional variant effect predictors that rely on predefined annotations, Evo2 learns patterns directly from genomic sequences.

The model uses a deep learning architecture called StripedHyena 2, which allows it to analyze long DNA sequences and learn relationships between nucleotides.

Evo2 can be used to investigate how genetic mutations affect sequence patterns learned during training, including mutations in non-coding regulatory regions.

This project uses GRCh38 (hg38) and the Evo2-40B model accessed through the NVIDIA API.

## How does Evo2 score variants?

Evo2 evaluates a variant by comparing the model's predicted likelihood of the reference nucleotide with that of the alternative nucleotide.

For each single-nucleotide variant (SNV), our scoring pipeline:

1. Extracts a 4,096-base-pair sequence surrounding the variant from the GRCh38 reference genome.  
2. Confirms that the reference nucleotide matches the reference genome.  
3. Uses Evo2 to obtain logits for the reference and alternative nucleotides.  
4. Calculates the difference between alternative and reference logits.  
5. Repeats the calculation using the reverse-complement sequence.  
6. Averages the forward and reverse scores to obtain the final Evo2 score.

The scoring formula is:

**Forward score:**

`Evo2_Forward = ALT_Logit − REF_Logit`

**Reverse score:**

`Evo2_Reverse = RC_ALT_Logit − RC_REF_Logit`

**Final score:**

`Evo2_Mean = (Evo2_Forward + Evo2_Reverse) / 2`

A negative score indicates that the alternative nucleotide is less favored by the model than the reference nucleotide in the evaluated sequence context. A positive score indicates that the alternative nucleotide is more favored.

## What is the output of Evo2?

Our Evo2 pipeline produces a comma-separated file (`.csv`) containing the variant information, experimental MPRA values, and Evo2 prediction scores.

The primary scoring file is:

`evo2_40b_api_checkpoint.csv`

| Column | Meaning |
| :---- | :---- |
| Variant\_ID | Unique genomic variant identifier |
| Chromosome | Chromosome containing the variant |
| Position | 1-based genomic position in GRCh38 |
| Ref | Reference nucleotide |
| Alt | Alternative nucleotide |
| Element | Regulatory element associated with the variant |
| MPRA\_Value | Experimental MPRA effect measurement |
| P\_Value | Statistical significance of the MPRA measurement |
| GRCh38\_REF | Reference nucleotide verified against GRCh38 |
| Evo2\_Forward\_REF\_Logit | Reference nucleotide logit in the forward orientation |
| Evo2\_Forward\_ALT\_Logit | Alternative nucleotide logit in the forward orientation |
| Evo2\_Forward | Forward ALT minus REF logit difference |
| Evo2\_Reverse\_REF\_Logit | Reference nucleotide logit in the reverse-complement orientation |
| Evo2\_Reverse\_ALT\_Logit | Alternative nucleotide logit in the reverse-complement orientation |
| Evo2\_Reverse | Reverse-complement ALT minus REF logit difference |
| Evo2\_Mean | Average of forward and reverse scores |
| Timestamp | UTC time when the variant was scored |

### **Notes**

* Higher positive values mean the alternative nucleotide is more favored relative to the reference; lower negative values mean it is less favored.  
* Evo2 scores are not PHRED-scaled and do not have established pathogenicity thresholds.  
* The same genomic SNV can occur in multiple MPRA elements. The sequence-based Evo2 score is calculated once per unique genomic variant and can be matched back to multiple experimental rows.  
* Failed scoring attempts are recorded separately in `evo2_40b_api_errors.csv`.  
* Successfully scored variants are saved continuously so interrupted jobs can resume without repeating completed predictions.

## How to run Evo2 with your own dataset:

### 1. Prepare your variants

Confirm that the dataset uses the GRCh38 (hg38) reference genome.

The input dataset should include:

* `Chromosome`  
* `Position`  
* `Ref`  
* `Alt`

The `Element`, `Value`, and `P-Value` columns are retained when available.

Only SNVs with single A, C, G, or T reference and alternative nucleotides are eligible for this pipeline.

Our project contains **29,444 SNV rows**, extracted from the larger dataset of 44,647 experimental variant rows.

The input file is:

`results/evo2_stratified_snvs.csv`

### 2. Prepare the reference genome

Download or provide the GRCh38 reference FASTA and ensure that it is accessible from the computing environment.

Our pipeline uses:

`reference/hg38.fa`

For every variant, the script extracts 2,048 bases upstream, the reference nucleotide, and 2,047 bases downstream.

The variant is evaluated using its surrounding genomic sequence.

Reference nucleotide mismatches and ambiguous sequence contexts are recorded as errors rather than silently scored.

### 3. Configure Evo2 through the NVIDIA API

Our project uses the NVIDIA-hosted Evo2 API rather than loading the entire model directly onto Kodiak GPUs.

The scoring script is:

`evo2_api_kodiak.py`

A valid NVIDIA API key must be available through the environment variable:

`NVIDIA_API_KEY`

On Kodiak, the key is loaded from a private environment file rather than being written directly into the Python script.

The script sends reference sequence contexts to the API and retrieves nucleotide logits used to calculate forward and reverse-complement scores.

API access, rate limits, and usage quotas should be checked before submitting large scoring jobs.

### 4. Submit the scoring job on Kodiak

Our scoring pipeline uses the PBS job scheduler.

The job script is:

`evo2_api_kodiak.pbs`

Submit the job using:

```
qsub evo2_api_kodiak.pbs
```

Check the job status:

```
qstat -u "$USER"
```

The PBS script runs the scoring program with the following arguments:

```
python evo2_api_kodiak.py \
  --input results/evo2_stratified_snvs.csv \
  --fasta reference/hg38.fa \
  --output results/evo2_40b_api_checkpoint.csv \
  --errors results/evo2_40b_api_errors.csv \
  --limit 0
```

The `--limit 0` option instructs the program to attempt all remaining eligible variants.

The current PBS configuration requests two CPU cores, 8 GB of memory, and an eight-hour walltime. Longer scoring runs may require multiple submissions.

### 5. Monitor scoring progress

The checkpoint file is updated after each successfully scored variant.

To count completed variants:

```
wc -l results/evo2_40b_api_checkpoint.csv
```

Subtract one for the header row.

To calculate progress:

```
import pandas as pd

df = pd.read_csv("results/evo2_40b_api_checkpoint.csv")

completed = df["Variant_ID"].nunique()
total = 29444

print(f"Completed: {completed:,} / {total:,}")
print(f"Remaining: {total-completed:,}")
print(f"Progress: {completed/total*100:.2f}%")
```

If the PBS job ends before scoring is complete, it can be resubmitted using the same checkpoint file. The script skips previously completed genomic variants.

Failed variants should be reviewed separately because they are not automatically considered completed.

### 6. Organize final results

Once scoring is complete, results will be organized in a folder named `RESULTS`.

Planned outputs include:

* `Evo2_All_Scores.tsv` — individual Evo2 predictions and variant metadata.  
* `FINAL_Combined_All_Models.tsv` — Evo2, AlphaGenome, CADD, AVI, MPRA, and conservation data.  
* `Evo2_Conservation/` — separate files for conservation groups such as mammals, primates, great apes, apes, human-chimp, and human-only.  
* `Evo2_Scoring_Summary.txt` — completed variants, missing scores, and scoring errors.

These outputs will support comparisons of Evo2 predictions against experimental MPRA effects and conservation levels.

&nbsp;