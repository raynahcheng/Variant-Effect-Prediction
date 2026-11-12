## About CADD

CADD (Combined Annotation Dependent Depletion) predicts how likely a human genetic variant is to be harmful or functionally disruptive. It is a machine learning model trained to tell apart two groups of variants, proxy-neutral and proxy-deleterious:

- Proxy-neutral: differences that became fixed on the human lineage since the human-chimp ancestor, so they are largely tolerated.
- Proxy-deleterious: simulated variants that mostly haven't been through selection.

This project uses GRCh38 (hg38), CADD v1.7.

## How does CADD score variants?

CADD combines annotations into one number per variant, including:

- Conservation (phyloP, GERP)
- Regulatory and epigenomic marks
- Gene model and consequence information
- Predicted protein effects

## What is the output of CADD?

CADD returns a gzip-compressed, tab-separated file (`.tsv.gz`). The first line is a comment (`##CADD GRCh38-v1.7`) and the second is the column header.

| Column | Meaning |
|---|---|
| Chrom, Pos, Ref, Alt | The variant, in VCF representation |
| RawScore | Direct model output. Higher means more deleterious. Only comparable within the same CADD version. |
| PHRED | Rank-based score. 10 = top 10%, 20 = top 1%, 30 = top 0.1%. Better for reporting and plotting. |


## How to run CADD with your own dataset:

### 1. Prepare your variants

Confirm your genome build (GRCh38 here). A wrong build makes every score wrong. Remove duplicate variants, since CADD only needs each one once.

### 2. Convert to VCF

Write a VCF (v4.2) with these rules:

- Columns: `#CHROM POS ID REF ALT QUAL FILTER INFO`, with `.` for ID, QUAL, FILTER and INFO.
- Chromosome names have no `chr` prefix (`1`, `2`, `X`).
- Sort by chromosome, then position.
- **Deletions need an anchor base.** `Alt = "-"` is not valid VCF. Use the base just before the deleted base as the anchor, and set POS to the anchor. Example: deleting C at 1000, where 999 is A, becomes `POS 999, REF AC, ALT A`.
- Best practice is to take anchor bases from the reference FASTA and left-normalize with `bcftools norm`.

### 3. Run CADD and match results

1. Go to https://cadd.gs.washington.edu/score and select GRCh38-v1.7.
2. Upload your VCF, choose whether to include annotations, and submit.
3. Download the `.tsv.gz` when the job finishes.
4. De-duplicate rows, then match to your original data on chromosome, position, ref and alt, applying the same deletion anchoring.
5. Check for variants with no score. Some edge-of-element deletions may be dropped.

For very large sets or fully reproducible runs, CADD can be installed locally (see the CADD-scripts repo on GitHub). Web queue times and file-size limits change, so check the site.
