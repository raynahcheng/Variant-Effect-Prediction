## About AlphaGenome Atlas

AlphaGenome Atlas is a precomputed resource from Google DeepMind containing predicted effects for every possible single-nucleotide variant in the human genome (about 9 billion). Instead of running a model yourself, you look up scores that have already been computed. For each variant the Atlas holds:

- Thousands of molecular effect predictions from AlphaGenome, covering gene regulation across hundreds of human and mouse cell types and tissues
- The AlphaGenome Variant Impact (AVI) score, a single number summarizing the variant's predicted impact
- Feature attributions that break the AVI score down into the biological features driving it
- Locations of over 2,500 DNA sequence motifs

The AVI score applies to both coding and non-coding regions. This project uses GRCh38 (hg38) coordinates.

## How does AlphaGenome Atlas score variants?

AVI combines several sources of evidence into one score per variant:

- AlphaGenome predictions of how the variant changes molecular readouts (expression, chromatin accessibility, splicing, TF binding and more)
- AlphaMissense, DeepMind's model for protein-altering (missense) variants
- Conservation across species (e.g. the Cactus alignment)

The score is decomposed into additive contributions from 18 modalities (e.g. Splicing, AlphaMissense, ChIP-TF, DNASE-seq, Cactus), so you can see which kind of evidence drives each variant's score. The top contributing modality is reported for every variant.

- The Atlas only covers single-nucleotide variants. Deletions are not scored, so in our table they have `AVI_SCORE` = 0 and are excluded from the AlphaGenome Atlas results on the Model Analytics page.
- Scores are precomputed, so there are no window-size or scorer choices to make, unlike running AlphaGenome directly.

## What is the output of AlphaGenome Atlas?

Each variant comes back with one row (or one VCF record) containing:

| Column | Meaning |
| --- | --- |
| `variant`, `chromosome`, `position`, `ref`, `alt` | The variant (1-based position) |
| `avi_raw` | The raw, uncalibrated model output. Higher means more impact. This is the `AVI_SCORE` column in our dataset. |
| `avi_phred` | Calibrated, PHRED-scaled score from 0 to about 70. 10 = top 10%, 20 = top 1%, 30 = top 0.1% of all genome-wide SNVs. Better for reporting. |
| `avi_quantile` | Calibrated tail quantile (1 − CDF), between 0 and 1 |
| `top_percentile` | Share of genome-wide SNVs with an equal or higher predicted impact |
| `top_modality`, `top_feature_importance` | The modality contributing most to the score, and its contribution |
| `fi_<MODALITY>` | Optional: the contribution of each of the 18 modalities |

In VCF output the same values are written as INFO tags and in a CSQ field (`Allele|AVI_PHRED|AVI_RAW|AVI_QUANTILE|AVI_TOP_PERCENTILE|AVI_TOP_FEATURE`).

## How to run AlphaGenome Atlas with your own dataset:

 ### &#32;1. Get an API key

The Atlas uses the same AlphaGenome API key as AlphaGenome (free for non-commercial use). Set it in your environment:

```bash
export ALPHAGENOME_API_KEY=your_key_here
```

For a handful of variants, you can skip the API and look them up on the [AlphaGenome Atlas website](https://alphagenome.google/atlas).

 ### 2. Prepare your variants

Use the same tab-separated file as for AlphaGenome:

| Column | Example | Notes |
| --- | --- | --- |
| `Chromosome` | `1` or `chr1` | Add the `chr` prefix if needed |
| `Position` | `109274968` | 1-based, GRCh38 |
| `Ref` | `G` | Reference allele |
| `Alt` | `T` | Alternate allele. Deletions (`-`) are not scored. |

Single variants are written as `chr:pos:ref>alt`, e.g. `chr1:109274968:G>T`.

 ### 3. Query the Atlas

Request the `AVI_SCORE` scorer (and `AVI_SCORE_FEATURE_IMPORTANCE` if you want the per-modality breakdown) for each variant through the AlphaGenome API. Batch lookups run at about 10 variants per second, so write large results straight to a file rather than holding them in memory.

 ### 4. Match results

1. Combine the results into one table with one row per variant.
2. Join back to your original data on `Chromosome`, `Position`, `Ref` and `Alt`.
3. Decide how to handle deletions and any other unscored variants. Leave them missing rather than filling with 0, so they don't count as low-impact variants.

For more detail, see the [AlphaGenome Atlas announcement](https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/), the [paper](https://storage.googleapis.com/deepmind-media/DeepMind.com/Blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/alphagenome-atlas.pdf) and the [AlphaGenome documentation](https://www.alphagenomedocs.com/).&#32;
