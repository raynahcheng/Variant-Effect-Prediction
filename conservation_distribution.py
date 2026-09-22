import pandas as pd

DATA_DIR = "/data/bentonm_shared/variant_effect_prediction/Variant-Effect-Prediction/data"

ORDER = ["human_only", "human_chimp", "great_apes", "apes", "primates", "mammals"]

df = pd.read_csv(f"{DATA_DIR}/variants_scored_conservation.tsv", sep="\t")
df["conservation_level"] = pd.Categorical(
    df["conservation_level"], categories=ORDER, ordered=True
)

print(df["conservation_level"].value_counts().sort_index())

for level in ORDER:
    subset = df[df["conservation_level"] == level]
    subset.to_csv(f"{DATA_DIR}/conservation_{level}.csv", index=False)
    print(f"Saved {len(subset)} variants to conservation_{level}.csv")