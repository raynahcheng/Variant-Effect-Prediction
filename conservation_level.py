import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor

DATA_DIR = "/data/bentonm_shared/variant_effect_prediction/Variant-Effect-Prediction/data"

# ── Species to check, mapped to their name in the cactus447way alignment ────
# (scientific-name prefix used in the MAF "src" field, e.g. "Pan_troglodytes.CM009238.2")
# No Nomascus_leucogenys (nomLeu3) in the 447-way tree, so gibbon is
# represented by Symphalangus_syndactylus (siamang) instead.
SPECIES_HIERARCHY = {
    # Great apes
    "chimp":      "Pan_troglodytes",
    "bonobo":     "Pan_paniscus",
    "gorilla":    "Gorilla_gorilla",
    "orangutan":  "Pongo_abelii",
    # Lesser apes
    "gibbon":     "Symphalangus_syndactylus",
    # Old World monkeys
    "macaque":    "Macaca_mulatta",
    "baboon":     "Papio_anubis",
    # New World monkeys
    "marmoset":   "Callithrix_jacchus",
    # Strepsirrhines
    "bushbaby":   "Otolemur_garnettii",
    # Non-primate mammals
    "mouse":      "Mus_musculus",
    "rat":        "Rattus_norvegicus",
    "dog":        "Canis_lupus_familiaris",
    "cow":        "Bos_taurus",
    "elephant":   "Loxodonta_africana",
}

CLADES = [
    ("human_chimp", ["chimp", "bonobo"]),
    ("great_apes",  ["chimp", "bonobo", "gorilla", "orangutan"]),
    ("apes",        ["chimp", "bonobo", "gorilla", "orangutan", "gibbon"]),
    ("primates",    ["chimp", "bonobo", "gorilla", "orangutan", "gibbon",
                     "macaque", "baboon", "marmoset", "bushbaby"]),
    ("mammals",     list(SPECIES_HIERARCHY.keys())),
]

# ── Query the cactus447way multiple alignment via UCSC's REST API ──────────
# Same alignment the phyloP447 scores (phyloP_scoring.py) are computed from,
# so "conserved at the species level" and "high phyloP" are consistent by
# construction instead of being two independently-noisy measurements
# (the old approach ran 14 separate pairwise hg38→species liftOvers via
# chain files, which can lose coverage at a locus for reasons unrelated to
# conservation: indels breaking a chain, repeat-masked regions, assembly
# gaps in the target genome, or short blocks dropped during chain filtering).
API_URL = "https://api.genome.ucsc.edu/getData/track"
GENOME = "hg38"
TRACK = "cactus447way"

session = requests.Session()

def species_aligned_at(chrom, pos_1based):
    """Returns the set of SPECIES_HIERARCHY keys with an aligned, non-gap base
    at this hg38 position in the cactus447way alignment."""
    chrom = str(chrom)
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    pos_0based = pos_1based - 1

    params = {
        "genome": GENOME, "track": TRACK,
        "chrom": chrom, "start": pos_0based, "end": pos_1based,
    }
    for attempt in range(3):
        try:
            resp = session.get(API_URL, params=params, timeout=30)
            resp.raise_for_status()
            blocks = resp.json().get(TRACK, [])
            break
        except requests.RequestException:
            if attempt == 2:
                raise

    prefix_to_species = {v: k for k, v in SPECIES_HIERARCHY.items()}
    present = set()

    for block in blocks:
        if not (block["chromStart"] <= pos_0based < block["chromEnd"]):
            continue
        offset = pos_0based - block["chromStart"]

        # mafBlock is a MAF alignment with newlines flattened to ';'. The
        # reference (hg38) row has no gaps in these reference-projected
        # blocks, so `offset` indexes directly into every row's text.
        for line in block["mafBlock"].split(";"):
            line = line.strip()
            if not line.startswith("s "):
                continue
            fields = line.split()
            src, text = fields[1], fields[-1]
            name = src.split(".")[0]
            if name == "hg38" or name not in prefix_to_species:
                continue
            if offset < len(text) and text[offset] not in "-Nn":
                present.add(prefix_to_species[name])
        break

    return present

def assign_conservation_level(coverage):
    level = "human_only"
    for clade_name, members in CLADES:
        if all(sp in coverage for sp in members):
            level = clade_name
    return level

df = pd.read_csv(f"{DATA_DIR}/variants_scored.tsv", sep="\t")

positions = list({(row.Chromosome, row.Position) for row in df.itertuples()})
print(f"Querying cactus447way alignment for {len(positions)} unique positions "
      f"across {len(df)} variants...")

with ThreadPoolExecutor(max_workers=10) as pool:
    coverages = pool.map(lambda cp: species_aligned_at(*cp), positions)
    levels = {key: assign_conservation_level(cov) for key, cov in zip(positions, coverages)}

df["conservation_level"] = [
    levels[(row.Chromosome, row.Position)] for row in df.itertuples()
]

ORDER = ["human_only", "human_chimp", "great_apes", "apes", "primates", "mammals"]
df["conservation_level"] = pd.Categorical(
    df["conservation_level"], categories=ORDER, ordered=True
)

print("\nConservation level counts:")
print(df["conservation_level"].value_counts().sort_index())

df.to_csv(f"{DATA_DIR}/variants_scored_conservation.tsv", sep="\t", index=False)
print("Saved.")
