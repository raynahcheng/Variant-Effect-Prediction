## About Evo2

Evo2 is a DNA sequence model that learns patterns from genomic sequence. For our project, we use it to score variants from an MPRA dataset and ask whether those scores agree with experimentally measured regulatory effects. We compare Evo2 with CADD and AlphaGenome, both overall and across groups defined by evolutionary conservation.

An important distinction: **Evo2 scores reflect the model's preference for one DNA sequence over another.** They are not direct predictions of gene expression, disease risk, or pathogenicity.

## Data used in our analysis

Our starting point is the project's MPRA dataset, including variants with both significant and nonsignificant experimental effects. Keeping both groups lets us evaluate whether Evo2 can distinguish variants with statistically significant MPRA results, rather than testing only on variants already known to have an effect.

The dashboard organizes variants into four conservation groups: **great apes, apes, primates, and mammals**. These groups come from cross-species sequence alignments. A position counts as present in a species if it has an aligned base; the base does not need to match the human reference.

### Input fields

| Field | Purpose |
|:--|:--|
| Chromosome and position | Locate the variant in GRCh38 |
| REF and ALT | Identify the reference and alternate alleles |
| MPRA measurement | Experimental regulatory effect for comparison |
| MPRA p-value | Label significant effects for classification analysis |
| Conservation group | Filter results by evolutionary conservation |

The scoring approach described here is for **single-nucleotide variants (SNVs)**. Indels need separate handling and should not be assumed to work with the same scoring procedure. Column names in the actual data files may differ from the labels above.

## How we score variants

Before scoring, we retrieve the surrounding sequence from the **GRCh38 human reference genome** and check that the REF allele matches the reference at the reported position. This catches coordinate and allele mismatches before they enter the analysis.

For each eligible SNV, we calculate two scores:

1. **Forward:** the alternate-allele logit minus the reference-allele logit, using the upstream sequence context.
2. **Reverse complement:** the corresponding alternate-minus-reference logit difference, using reverse-complemented downstream context.

We take the average as the final score:

$$
\text{Evo2 mean score} = \frac{\text{forward score} + \text{reverse-complement score}}{2}
$$

A positive score means Evo2 favors the alternate allele under this scoring method; a negative score means it favors the reference allele. A larger absolute score indicates a stronger preference, **not necessarily a larger change in expression**. We keep the forward and reverse-complement scores separately so we can check how closely the two orientations agree.

## Output and matching

The Evo2 output retains the original variant identifiers and adds the forward score, reverse-complement score, and mean score. For analysis of effect magnitude, we can also calculate the absolute mean score.

We match Evo2 predictions to the shared dataset using the full variant key: **chromosome, position, REF, and ALT**. This matters because the same genomic variant can appear in more than one experimental element or cell type, with a different MPRA measurement for each record. We preserve those measurements rather than collapsing them into a single value.

Before plotting results, we check the number of matched and unmatched records and make sure the model comparisons use the same eligible variants wherever possible.

## Reading the dashboard plots

### Kendall's tau

Kendall's tau tells us how closely the rankings of Evo2 scores and MPRA measurements agree. For the magnitude-based comparison, we use the absolute mean Evo2 score and absolute MPRA effect. A higher positive tau means the two rankings are more consistent. It does not mean Evo2 predicts the direction of expression change.

### ROC curves and AUROC

We label an MPRA result significant when its **p-value is below 0.05**. We then use the absolute mean Evo2 score to rank variants and calculate the ROC curve and area under the curve (AUROC). An AUROC near 0.5 suggests little ability to distinguish the two classes; higher values indicate better separation on the evaluated dataset.

### Precision–recall

Precision–recall analysis is useful when the number of significant variants differs from the number of nonsignificant variants. We interpret AUPRC alongside the proportion of significant variants in each group.

### Conservation groups and model comparisons

The dashboard lets users compare models and select conservation groups. We interpret differences carefully because the groups do not all contain the same number or mix of variants. The same caution applies when comparing transposable and nontransposable elements.

## Using Evo2 with another dataset

This section outlines the workflow we used; it is **not a copy-and-run script**. Exact commands depend on the final scoring code and computing environment.

1. Prepare an SNV table with GRCh38 chromosome, position, REF, and ALT fields.
2. Retrieve GRCh38 sequence context and confirm that each REF allele matches the reference.
3. Run Evo2 scoring in the forward and reverse-complement orientations.
4. Calculate the mean score and save all three score columns alongside the original variant identifiers.
5. Join the scores back to experimental measurements using chromosome, position, REF, and ALT.
6. If MPRA measurements are available, calculate Kendall's tau and classification metrics using consistent definitions.
7. Add conservation-group annotations if you want to reproduce the dashboard's stratified comparisons.

For this project, we first tested Evo2 in Google Colab, then moved the workflow to Baylor’s Kodiak computing environment to support longer scoring runs on the full MPRA dataset. The same scoring and matching steps are used to prepare Evo2 results for comparison with CADD and AlphaGenome in our dashboard. Anyone adapting the workflow to another dataset should use the final validated scoring scripts and confirm the model version, dependencies, and sequence-context settings before running it.

## How Evo2 fits into our project

Our project asks whether variant-effect prediction models agree more closely with experimental MPRA measurements at genomic positions conserved across broader evolutionary groups. Evo2 contributes a sequence-based perspective alongside CADD and AlphaGenome. By joining all three sets of predictions to the same variants, we can compare their performance using Kendall’s tau, ROC curves, and conservation-group filters in the dashboard.

These comparisons also help us identify where model predictions differ from experimental results. Evo2 learns patterns in DNA sequence, whereas MPRA measures regulatory activity directly, so agreement is not guaranteed. We retain the forward and reverse-complement scores to examine orientation differences, and we interpret comparisons in light of the different sample sizes and variant compositions across conservation groups.

Our Evo2 workflow currently focuses on SNVs. As we finalize the merged dataset and plots, we will verify that the scored variants match the dashboard data before drawing conclusions about Evo2’s performance relative to the other models. This keeps the dashboard focused on our central research question: **how does evolutionary conservation relate to the accuracy of computational variant-effect predictions?**
