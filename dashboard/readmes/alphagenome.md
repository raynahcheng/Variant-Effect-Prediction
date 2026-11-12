## About AlphaGenome

AlphaGenome is a deep learning model from Google DeepMind that reads DNA sequence and predicts what that sequence does at the molecular level. Given up to about 1 Mb of input sequence, it predicts thousands of functional genomic tracks across many human (and mouse) cell types and tissues, including:

- Gene expression: RNA-seq, CAGE, PRO-cap
- Chromatin accessibility: DNase-seq, ATAC-seq
- Histone modifications and transcription factor binding: ChIP-seq
- Splicing: splice sites, splice site usage, splice junctions
- 3D genome structure: contact maps

Most tracks are predicted at single-base resolution. Histone and TF ChIP tracks are predicted at 128 bp resolution, and contact maps at 2,048 bp.

This project uses GRCh38 (hg38) coordinates and scores variants through the AlphaGenome API.

## How does AlphaGenome score variants?

AlphaGenome scores a variant by running the model twice on the same stretch of genome: once with the reference allele and once with the alternate allele. A variant scorer then compares the two predictions and summarizes the difference as a number.

- The variant sits at the center of a window of fixed length (2 KB, 16 KB, 100 KB, 500 KB or 1 MB). We use 1 MB, the largest window, so distal regulatory effects on nearby genes can be captured.
- Predictions are computed on demand by the API for each request.
- AlphaGenome provides a set of recommended variant scorers, one per output type. We use the RNA\_SEQ scorer because our questions is about effects on expression. This scorer estimates how the alternate allele changes predicted expression of each gene in the window, as a log fold change over that gene's exons.
- One variant therefore produces many scores: one per (gene × track), where tracks are the individual cell types and tissues.

## What is the output of AlphaGenome?

The script calls `variant_scorers.tidy_scores()` on the API results, which produces one row per **variant × gene × track**. Results are saved in chunks as Parquet files to reduce the amount of variants being called by the API at a time (`scores_chunk_00000.parquet`, `scores_chunk_00001.parquet`, …).

Key columns:

| Column | Meaning |
| --- | --- |
| `Chromosome`, `Position`, `Ref`, `Alt` | The variant as it appears in the input table, so rows can be joined back to the original data |
| `variant_id` | The variant as AlphaGenome represents it (VCF-style alleles, `chr` prefix) |
| `scored_interval` | The genomic window the model was run on |
| `gene_id`, `gene_name`, `gene_type`, `gene_strand` | The gene the score refers to |
| `output_type`, `variant_scorer` | Which output (e.g. RNA\_SEQ) and scorer produced the score |
| `track_name`, `biosample_name`, `ontology_curie`, … | Which cell type, tissue or experiment the track represents |
| `raw_score` | The direct scorer output (log fold change, alternate vs. reference). Sign gives direction. |
| `quantile_score` | The raw score ranked against a background distribution of variants for that track, putting tracks on a common scale. Better for comparing across tracks. |

## How to run AlphaGenome with your own dataset:

 ### &#32;1. Get an API key and install the package

Request an AlphaGenome API key from Google DeepMind (it is free for non-commercial use, and usage is rate-limited). Then install the Python client:

```bash
pip install alphagenome
export ALPHAGENOME_API_KEY=your_key_here
```

 ### 2. Prepare your variants

Create a tab-separated file with these columns:

| Column | Example | Notes |
| --- | --- | --- |
| `Chromosome` | `1` or `chr1` | The script adds the `chr` prefix if it's missing |
| `Position` | `109274968` | 1-based, GRCh38 |
| `Ref` | `G` | Reference allele |
| `Alt` | `T` or `-` | Alternate allele; `-` means a deletion |

 ### 3. Choose your scorer(s)

The default in our scripts is RNA\_SEQ. To use other track types, pass them with `--scorers` (any of AlphaGenome's recommended scorers, e.g. `DNASE`, `ATAC`, `CAGE`, `CHIP_HISTONE`). Pick the output type that best matches the biology you're measuring.

 ### &#32;4. Run the scoring script

Use [alphagenome\_scoring.py](https://github.com/raynahcheng/Variant-Effect-Prediction/blob/raynah/alphagenome_scoring.py):

```bash
python alphagenome_scoring.py --limit 5                  # quick test on 5 variants
python alphagenome_scoring.py                            # all variants
python alphagenome_scoring.py --skip-deletions           # SNVs only
python alphagenome_scoring.py --scorers DNASE CAGE       # other track types
python alphagenome_scoring.py --input my_variants.tsv --outdir results/
```

Useful options:

| Option | Default | Meaning |
| --- | --- | --- |
| `--input` | `data/variants_all_te.tsv` | Input variant table |
| `--outdir` | `data/alphagenome_scores_<scorers>` | Where Parquet chunks are written |
| `--sequence-length` | `1MB` | Window size around each variant |
| `--chunk-size` | `100` | Variants per output file |
| `--max-workers` | `5` | Concurrent API requests |

The run is restartable: finished chunks are skipped, so if the job crashes or hits a walltime limit, run the same command again. If one variant in a chunk fails, the script retries that chunk's variants one at a time so the rest still get scored.

### 5. Combine and match results

1. Read all chunks into one table, e.g. `pd.concat(pd.read_parquet(f) for f in glob.glob("outdir/scores_chunk_*.parquet"))`.
2. Aggregate to one row per variant (see Notes above).
3. Join back to your original data on `Chromosome`, `Position`, `Ref` and `Alt`.

For more detail on the API, available tracks and scorers, see the [AlphaGenome documentation](https://www.alphagenomedocs.com/colabs/v1/quick_start.html).&#32;
