# AlphaGenome

AlphaGenome is a Google DeepMind model that reads DNA sequence and predicts what that DNA does biochemically. It takes up to about a million bases of input and predicts many data types at single-base resolution. These include gene expression (RNA-seq, CAGE), chromatin accessibility (DNase, ATAC), histone marks, transcription factor binding, splicing and 3D contact maps. Each is predicted across many cell types and tissues.

For variant work, it predicts the same region twice, once with the reference base and once with your alternate base, and reports the difference. You reach it through the alphagenome Python package, which calls Google’s hosted service with a personal API key

## Overview
- **Input:** 
A CSV of variants on the GRCh38 genome build, with chromosome, position, ref and alt columns. The original file also carried three MPRA columns, the element name, the measured value and its p-value. The code adds a “chr” prefix to chromosome names and turns dashes (deletions) into empty strings. Each variant is scored in a 16 KB window centered on it, using four scorers: RNA-seq, CAGE, DNase and ATAC. The full file had 44,647 variants.


- **Output:** 
AlphaGenome returned its standard results file, which is very large because every variant is scored against many genes and cell types or tissues. We then used feature extraction to loop through that file and pull out the identifier and raw score for each unique variant. We joined those to the MPRA values and the other columns from the original file. The result is a much smaller table with one row per variant, holding the predicted score next to the experimental measurement, ready for statistics and plotting.


- **Visualizations:** The results are the Kendall tau statistics graph and the ROC curves, both built from the extracted AlphaGenome raw scores and the MPRA values. The Kendall tau graph shows how well the model’s predicted ranking of variants agrees with the ranking from the lab measurements. A tau near 1 means the rankings agree, a tau near 0 means there is no relationship, and a negative tau means they run opposite. The ROC curves show how well the predicted scores separate variants that had a real effect in the MPRA from those that did not, with the area under each curve summarizing that ability. I don’t have your actual tau values or AUC numbers, so add those here to complete this section


- **Why Only RNASeq scores where taken:** We used only RNA-seq scores because an MPRA measures how much a variant changes expression, and RNA-seq is the AlphaGenome output that predicts expression most directly. Using a single score type also kept each variant to one predicted value, which made the comparison cleaner.


- **Visualizations Decision Making:**
  - Kendall tau was chosen because it compares rankings instead of raw values. AlphaGenome raw scores are very small numbers on a different scale from MPRA values, so a measure that depends on ranks is more meaningful than one that expects the two to match numerically. It also holds up better against outliers and nonlinear relationships. The graph makes it easy to see at a glance whether the model orders variants the way the experiment did.

  - ROC curves were chosen to test whether the model can act as a classifier. They turn the question into whether the predicted score can tell apart variants with significant MPRA effects from those without. This complements tau, which looks at the whole ranking, by focusing on the practical question of picking out the variants that matter. One thing to confirm is how you defined a variant with a real effect, such as an MPRA p-value cutoff, since that definition sets the labels the ROC curve depends on, and it’s worth stating in this section.
