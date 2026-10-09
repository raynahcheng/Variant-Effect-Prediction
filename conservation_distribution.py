import argparse
import pandas as pd

DATA_DIR = "/data/bentonm_shared/variant_effect_prediction/Variant-Effect-Prediction/data"

ORDER = ["human_only", "human_chimp", "great_apes", "apes", "primates", "mammals"]

parser = argparse.ArgumentParser(description="Split variants into one CSV per conservation level.")
parser.add_argument("--input", default=f"{DATA_DIR}/variants_scored_conservation.tsv")
parser.add_argument("--prefix", default="conservation",
                    help="writes data/<prefix>_<level>.csv")
args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t")
df["conservation_level"] = pd.Categorical(
    df["conservation_level"], categories=ORDER, ordered=True
)

print(df["conservation_level"].value_counts().sort_index())

for level in ORDER:
    subset = df[df["conservation_level"] == level]
    subset.to_csv(f"{DATA_DIR}/{args.prefix}_{level}.csv", index=False)
    print(f"Saved {len(subset)} variants to {args.prefix}_{level}.csv")
