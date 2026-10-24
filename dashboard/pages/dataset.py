import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, callback, dcc, html

from data import DATA_PATH, VARIANTS

dash.register_page(__name__, path="/dataset", name="Dataset", order=3)

LINK_STYLE = {"color": "#1f4e8c"}
KIRCHER_URL = "https://www.nature.com/articles/s41467-019-11526-w"

layout = html.Div(
    [
        html.H2("Dataset"),
        html.P([
            "Our dataset of variants come from the saturation mutagenesis massively parallel reporter assay "
            "(MPRA) of ",
            html.A("Kircher et al. (2019)", href=KIRCHER_URL, target="_blank", style=LINK_STYLE),
            ", which was also used in the CAGI5 regulatory variant "
            f"challenge. In this experiment, every possible single nucleotide substitution, along with "
            f"many 1 bp deletions, was introduced into 21 disease-associated regulatory elements, and "
            f"the effect of each variant on reporter gene expression was measured. Ten of the elements "
            f"are promoters (F9, FOXE1, GP1BA, HBB, HBG1, HNF4A, LDLR, MSMB, PKLR, and TERT) and eleven "
            f"are enhancers (BCL11A, IRF4, IRF6, two MYC enhancers tagged by rs6983267 and rs11986220, "
            f"RET, SORT1, TCF7L2, UC88, ZFAND3, and the ZRS limb enhancer of SHH). Some elements were "
            f"assayed more than once, so they appear under several names in the Element column: TERT "
            f"was tested in four cell lines (TERT-GAa, -GBM, -GSc, -HEK), PKLR at two time points "
            f"(24h and 48h), LDLR and SORT1 in replicate experiments (.2), SORT1 in reverse orientation "
            f"(SORT1-flip), and two ZRS haplotypes (ZRSh-13, ZRSh-13h2). Altogether the dataset holds "
            f"{len(VARIANTS):,} variants across {VARIANTS['Element'].nunique()} element experiments. "
            f"Each variant's Value is the log2 fold change in expression relative to the reference "
            f"sequence."],
            style={"textAlign": "center", "maxWidth": "900px", "margin": "0 auto 1rem"},
        ),
        html.P([
            "The downloadable file contains the original ",
            html.A("Kircher et al.", href=KIRCHER_URL, target="_blank", style=LINK_STYLE),
            " dataset (Chromosome, Position, "
            "Ref, Alt, Tags, DNA, RNA, Value, P-Value, and Element) along with the annotations we added: "
            "phyloP conservation scores (phyloP447), conservation level annotations (conservation_level), "
            "the species each site aligns to in the UCSC Cactus 447-way alignment (aligned_species), "
            "transposable element status from RepeatMasker (te_status, rep_class, rep_family, rep_name, "
            "and rep_divergence), and raw scores from AlphaGenome (raw_score) and CADD (CADD_raw_score "
            "and CADD_PHRED). Evo2 scores will be added once they are available."],
            style={"textAlign": "center", "maxWidth": "900px", "margin": "0 auto 1rem"},
        ),
        html.Div(
            dbc.Button("Download Dataset", id="ds-download-btn", color="secondary"),
            className="text-center",
        ),
        dcc.Download(id="ds-download"),
    ]
)


@callback(
    Output("ds-download", "data"),
    Input("ds-download-btn", "n_clicks"),
    prevent_initial_call=True,
)
def download(n):
    return dcc.send_file(DATA_PATH)
