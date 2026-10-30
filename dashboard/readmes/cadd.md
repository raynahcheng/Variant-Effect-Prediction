
## 1. Project summary

We have a set of single-nucleotide substitutions and single-base deletions across several regulatory elements, including BCL11A, SORT1, TERT, MYC, IRF4, LDLR, UC88 and ZRSh-13. Each variant has experimental measurements from a massively parallel reporter assay, along with conservation and repeat annotations.

The variants were originally split into four CSV files by conservation level, meaning the set of species in which the position aligns: apes, great apes, primates and mammals. The apes file has 132 rows, great apes has 123, primates has 4,330, and mammals has 40,062, for 44,647 rows in total.

The goal was to obtain a CADD (Combined Annotation Dependent Depletion) score for every variant, so that we can compare predicted deleteriousness with the experimental results and with conservation.

Many genomic variants appear under more than one element, for example SORT1, SORT1-flip and SORT1.2, or several cell-type versions of TERT. CADD only needs each variant once, so after removing duplicates we submitted 95 variants for apes, 123 for great apes, 2,851 for primates and 30,445 for mammals.

---

## 2. Workflow at a glance

First, each CSV was converted to a VCF file in GRCh38 coordinates, with deletions rewritten in proper VCF form. Second, each VCF was uploaded to the CADD web service and the scored output was downloaded. The CADD output can then be matched back to the original data using chromosome, position, reference allele and alternate allele.

---

## 3. Genome build

We used GRCh38 (hg38) and CADD version 1.7.

The build was inferred, not stated in the source CSVs. The evidence is that the MYC rs6983267 and rs11986220 elements sit at roughly chr8:127.4 Mb and 127.5 Mb in our data. That matches GRCh38, whereas in GRCh37 those loci are near chr8:128.4 Mb. The BCL11A element near chr2:60.49 Mb is also consistent with GRCh38. The CADD result files confirm the version they were run with, since their header line reads "CADD GRCh38-v1.7".

It is worth confirming the build with whoever produced the original CSVs. If it turns out to be GRCh37, every score is wrong and the VCFs must be rescored against the GRCh37 model, or lifted over first.

---

## 4. Step 1: Converting the CSVs to VCF

### 4.1 Input format

Each CSV has, among other columns, Chromosome, Position, Ref and Alt. Substitutions have a single base in both Ref and Alt. Deletions have Alt set to "-". There are no insertions in these data.

### 4.2 VCF format used

The VCFs are version 4.2. Chromosome names have no "chr" prefix (1, 2, X and so on), which is what CADD expects. The ID, QUAL, FILTER and INFO columns are all ".", so no metadata is carried through. Records are sorted by chromosome and then position, and exact duplicates are removed. A short example:

```
##fileformat=VCFv4.2
##reference=GRCh38 (assumed)
#CHROM  POS        ID  REF  ALT  QUAL  FILTER  INFO
1       109274976  .   T    A    .     .       .
1       109274976  .   TC   T    .     .       .
```

### 4.3 How deletions were handled

VCF cannot represent a deletion as Alt equal to "-". A deletion in VCF must include an anchor base, which is the base immediately before the deleted base, in both REF and ALT, and POS must point at the anchor.

For example, suppose the base C at position 1000 is deleted and the base at position 999 is A. The CSV describes this as position 1000, Ref C, Alt "-". The VCF describes it as position 999, REF "AC", ALT "A".

Where the anchor base came from: we did not have a reference genome available in the working environment, so we built a lookup of known reference bases from the Ref column of the input rows themselves, mapping each chromosome and position to its reference base. The element data is a dense tiling of positions, so the base at the position just before a deletion is almost always present as the Ref of some other row. We checked the lookup for conflicts, meaning the same position reported with two different reference bases, and found none.

Fallback for element edges: for deletions at the very first position of an element, the preceding base is not in the data. For these we anchored on the following base instead. To delete the base at position p, we write POS as p, REF as the base at p followed by the base at p+1, and ALT as the base at p+1. This is a valid VCF representation of the same deletion, but it is not left-aligned. CADD normalizes indels internally, so this should be fine, but it is the less preferred route.

Across all four files, before removing duplicates across elements, 5,452 deletions were anchored on the preceding base, 25 used the following-base fallback, and none could not be anchored.

One caveat: the anchor bases were never checked against the real reference genome. If the CSV Ref values were wrong, or the build is wrong, CADD may have rejected or mis-scored the affected deletions (see section 7).

### 4.4 Conversion script

This is `make_vcf.py`:

```python
import pandas as pd

UP = '/mnt/user-data/uploads'
OUT = '/mnt/user-data/outputs'
files = {
    'apes':       'conservation_all_apes.csv',
    'great_apes': 'conservation_all_great_apes.csv',
    'mammals':    'conservation_all_mammals.csv',
    'primates':   'conservation_all_primates.csv',
}
dfs = {k: pd.read_csv(f'{UP}/{v}', keep_default_na=False, dtype={'Chromosome': str})
       for k, v in files.items()}

# Reference-base lookup built from the Ref column of every row in all files
refmap = {}
for df in dfs.values():
    for c, p, r in zip(df.Chromosome, df.Position, df.Ref):
        assert refmap.get((c, p), r) == r, f'conflicting ref at {c}:{p}'
        refmap[(c, p)] = r

order = {str(i): i for i in range(1, 23)}
order.update({'X': 23, 'Y': 24})
header = ['##fileformat=VCFv4.2',
          '##reference=GRCh38 (assumed)',
          '##source=converted from MPRA conservation CSV',
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO']

for name, df in dfs.items():
    records = set()
    for c, p, r, a in zip(df.Chromosome, df.Position, df.Ref, df.Alt):
        if a != '-':                                   # substitution
            records.add((c, p, r, a))
        elif (c, p - 1) in refmap:                     # deletion, left anchor
            b = refmap[(c, p - 1)]
            records.add((c, p - 1, b + r, b))
        elif (c, p + 1) in refmap:                     # deletion, right-anchor fallback
            b = refmap[(c, p + 1)]
            records.add((c, p, r + b, b))
        else:
            raise ValueError(f'cannot anchor deletion at {c}:{p}')
    rows = sorted(records, key=lambda x: (order[x[0]], x[1], x[2], x[3]))
    with open(f'{OUT}/{name}.vcf', 'w') as f:
        f.write('\n'.join(header) + '\n')
        for c, p, r, a in rows:
            f.write(f'{c}\t{p}\t.\t{r}\t{a}\t.\t.\t.\n')
    print(name, len(df), '->', len(rows))
```

---

## 5. Step 2: Running CADD

### 5.1 What CADD is

CADD scores how deleterious, meaning how likely to be functionally disruptive, a human variant is predicted to be. It is a machine-learning model trained to tell apart two kinds of variants. The first are proxy-neutral variants, which are differences that arose and became fixed on the human lineage since the human-chimp ancestor and are therefore largely tolerated by selection. The second are proxy-deleterious variants, which are simulated variants that mostly have not been through selection.

CADD integrates dozens of annotations into a single number per variant. These include conservation scores such as phyloP and GERP, regulatory and epigenomic marks, gene model and consequence information, and predicted protein effects.

### 5.2 Running it through the web service

Go to the CADD scoring site at https://cadd.gs.washington.edu/score and select the genome build and version GRCh38-v1.7. Choose one of the four VCF files as the input, decide whether to include annotations, submit, wait for the job to finish, and download the resulting .tsv.gz file. Repeat for each VCF.

On annotations, our apes and great apes results contain only six columns (Chrom, Pos, Ref, Alt, RawScore, PHRED), so they are scores only. Our primates and mammals results contain the full annotation output of about 153 columns. Either choice works for this project, because we only use the raw score and the PHRED score.

Web-service queue times and file-size limits change over time, so check the site for current limits. For very large sets, or for a fully reproducible local run, CADD can also be installed and run locally (see the CADD-scripts repository on GitHub). The local version needs the large annotation download and the matching reference genome.

### 5.3 Output format

The result files are tab-separated and gzip-compressed. The first line is a comment beginning with "##CADD GRCh38-v1.7", and the second line is the column header, beginning with "#Chrom".

The Chrom, Pos, Ref and Alt columns are in VCF representation, so deletions appear as, for example, "TC to T" rather than "C to -". Any matching back to the original data therefore has to apply the same anchoring rule described in section 4.3.

Scores-only files have exactly one row per variant. Annotated files have several rows per variant, one for each overlapping transcript, consequence or feature, with identical RawScore and PHRED values on every row. These must be de-duplicated before matching. We verified that across all 33,510 distinct scored variants, no variant had conflicting scores between rows or between files.

---

## 6. How to interpret CADD scores

### 6.1 Two numbers per variant

The raw score is the direct output of the model. Higher means the variant looks more like the simulated, deleterious class, and lower, including negative values, means it looks more like the proxy-neutral class. The absolute value only has meaning relative to other variants scored by the same CADD version. Raw scores are not comparable across CADD versions or genome builds.

The PHRED-scaled score is a rank-based transformation of the raw score. It compares the variant to all possible single-nucleotide substitutions in the human genome, and is defined as minus 10 times the base-10 logarithm of the fraction of all possible substitutions that score at least as high. In practice, a PHRED of 10 means the variant is in the top 10 percent of possible substitutions, 20 means the top 1 percent, 30 means the top 0.1 percent, and 40 means the top 0.01 percent.

A common rule of thumb treats a PHRED of 10 or more as notable and 20 or more as likely highly impactful, and a cutoff of about 15 is often used for coding variants. These are heuristics and not validated thresholds for MPRA effect size. For comparisons within this project, rank-based reasoning is safer than a hard cutoff.

Because PHRED is a rank, it is easier to compare across variants and across CADD versions, so it is the better choice for reporting and plotting. The raw score gives finer resolution among low-scoring variants, since PHRED compresses everything below about 1 into values near zero.

### 6.2 What CADD does and does not tell you

CADD predicts general deleteriousness, which is closely related to evolutionary constraint. It was not trained on MPRA expression effects. A high CADD score therefore does not imply a large change in reporter expression, and a low score does not imply no change.

CADD's input features include conservation scores such as phyloP and GERP, so CADD scores will correlate with the phyloP and conservation-level information in this project. That correlation is partly circular. When asking whether conservation predicts regulatory effect, remember that CADD already contains conservation information and is not an independent line of evidence about conservation. It is better to compare the experimental effect against phyloP and against CADD separately than to treat CADD as a second, independent predictor.

For deletions, CADD uses an indel-specific scoring path. Indel scores are generally less well calibrated than substitution scores, so a single-base deletion's PHRED should not be compared one-to-one with a substitution's PHRED without caution.

Variants in repetitive or transposable-element sequence can be harder to score reliably because of alignment and mapping issues, which is relevant to our repeat annotations.

### 6.3 Suggested analyses

Reasonable next steps include correlating the CADD PHRED score with the experimental effect, split by conservation level and by element. It would also be useful to check whether substitutions and deletions behave differently, and to compare high-PHRED variants that show no experimental effect against low-PHRED variants that do.

---

## 7. Known issues and limitations

The genome build is inferred and not confirmed (section 3).

Four distinct deletions received no score from CADD. All are single-base deletions at the edge of an element, and all are in the mammals set. They are the deletion of G at chr1:109274652 (SORT1 and SORT1-flip), the deletion of T at chr1:109274651 (SORT1-flip), the deletion of A at chr2:161238998 (UC88), and the deletion of A at chr7:156791604 (ZRSh-13 and ZRSh-13h2). In the original tables these account for six rows. All four were present in the submitted VCFs, but CADD returned no output for them. The cause is unknown. One possibility is a mismatch between our inferred anchor base and the true reference. To investigate, extract the corresponding VCF lines, check the anchor and deleted bases against the real GRCh38 sequence, and resubmit.

The anchor bases were inferred from the data and not from the reference genome. A cleaner pipeline would left-normalize the deletions against a reference FASTA, for example with bcftools norm, before submitting.

The 25 right-anchored deletions are not left-aligned. CADD should normalize them, but it is worth verifying if any of them score oddly.

Scores are per variant and not per element. The same genomic variant tested in several elements gets the same CADD score each time. That is correct, since CADD only sees the genome, but it means those rows are not independent when running statistics.
