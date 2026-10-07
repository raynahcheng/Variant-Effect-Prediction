import argparse
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor

DATA_DIR = "/data/bentonm_shared/variant_effect_prediction/Variant-Effect-Prediction/data"

parser = argparse.ArgumentParser(description="Assign cactus447way conservation levels.")
parser.add_argument("--input", default=f"{DATA_DIR}/variants_all_phyloP.tsv")
parser.add_argument("--output", default=f"{DATA_DIR}/variants_all_conservation.tsv")
args = parser.parse_args()

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

# Each level lists only the species that branch off at that level (closest to
# most distant from human). A variant gets the most distant level at which AT
# LEAST ONE species has an aligned base. Requiring every species instead let a
# single gap or unsequenced (N) base in one genome demote a site, e.g. bonobo
# is N across the TERT promoter, which put sites aligned in chimp, gorilla,
# orangutan, macaque, etc. into human_only.
LINEAGES = [
    ("human_chimp", ["chimp", "bonobo"]),
    ("great_apes",  ["gorilla", "orangutan"]),
    ("apes",        ["gibbon"]),
    ("primates",    ["macaque", "baboon", "marmoset", "bushbaby"]),
    ("mammals",     ["mouse", "rat", "dog", "cow", "elephant"]),
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

        # mafBlock is a MAF alignment with newlines flattened to ';'.
        rows = {}
        for line in block["mafBlock"].split(";"):
            line = line.strip()
            if line.startswith("s "):
                fields = line.split()
                rows[fields[1].split(".")[0]] = fields[-1]

        # The hg38 row CAN contain gaps (columns where another species has an
        # insertion), so the alignment column for this base is the offset-th
        # non-gap character of the hg38 row, not simply `offset`.
        hg_cols = [i for i, c in enumerate(rows["hg38"]) if c != "-"]
        col = hg_cols[offset]

        # "-" = gap (no aligned base). N = base unknown in that assembly, which
        # says nothing about whether the sequence exists, so it isn't counted
        # as present either; under the at-least-one rule below it can no
        # longer demote a site on its own.
        for name, text in rows.items():
            if name in prefix_to_species and text[col] not in "-Nn":
                present.add(prefix_to_species[name])
        break

    return present

def assign_conservation_level(coverage):
    level = "human_only"
    for lineage_name, members in LINEAGES:
        if any(sp in coverage for sp in members):
            level = lineage_name
    return level

df = pd.read_csv(args.input, sep="\t")

positions = list({(row.Chromosome, row.Position) for row in df.itertuples()})
print(f"Querying cactus447way alignment for {len(positions)} unique positions "
      f"across {len(df)} variants...")

with ThreadPoolExecutor(max_workers=10) as pool:
    coverages = pool.map(lambda cp: species_aligned_at(*cp), positions)
    coverage_by_pos = dict(zip(positions, coverages))

df["conservation_level"] = [
    assign_conservation_level(coverage_by_pos[(row.Chromosome, row.Position)])
    for row in df.itertuples()
]
# Which species aligned, so a level can be checked by hand.
df["aligned_species"] = [
    ",".join(sp for sp in SPECIES_HIERARCHY if sp in coverage_by_pos[(row.Chromosome, row.Position)])
    for row in df.itertuples()
]

ORDER = ["human_only", "human_chimp", "great_apes", "apes", "primates", "mammals"]
df["conservation_level"] = pd.Categorical(
    df["conservation_level"], categories=ORDER, ordered=True
)

print("\nConservation level counts:")
print(df["conservation_level"].value_counts().sort_index())

df.to_csv(args.output, sep="\t", index=False)
print("Saved.")
