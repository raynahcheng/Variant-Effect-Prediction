## About AlphaGenome
AlphaGenome is a Google DeepMind deep learning model that reads DNA sequence and predicts what that DNA does biochemically. It takes up to about a 1 Mb of input and predicts many tracks at single-base resolution. These include gene expression (RNA-seq, CAGE), chromatin accessibility (DNase, ATAC), histone marks, transcription factor binding, splicing and 3D contact maps. Each is predicted across many cell types and tissues.

## How does AlphaGenome score variants?
For variant work, it predicts the same region twice, once with the reference base and once with your alternate base, and reports the difference. AlphaGenome has already calculated these tracks for each variant, so calling each variant with a personal API key retrieves the prediction for that variant. 

## What is the output of AlphaGenome?
Depending on which tracks you select, a variant will be scored based on those specific tracks with a raw score, a absolute maximum, and an absolute minimum score. We chose to utilize the RNA-seq track because this made the most sense for our data. AlphaGenome outputs a results file with however many variants you called with the API key, and its track scores.

## How to run AlphaGenome with your own dataset:
1. Create a CSV containing the variants you want to predict, including columns for chromosome, position in the chromosome, reference allele, and the alternative allele. 
2. Choose which track is best to run for your variants. If you choose more than one, it is best to aggregate the raw scores and take the average to get the best result for your model prediction.
3. Obtain the variant prediction scores using the AlphaGenome API key.
4. Utilize the results file to visualize variant effects.

The documentation for AlphaGenome is [linked here](https://www.alphagenomedocs.com/colabs/v1/quick_start.html) for your reference. To score our own variants, we used this script: [alphagenome_scoring.py](https://github.com/raynahcheng/Variant-Effect-Prediction/blob/raynah/alphagenome_scoring.py).

