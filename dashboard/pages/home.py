import dash
import dash_bootstrap_components as dbc
from dash import dcc, html

dash.register_page(__name__, path="/", name="Home", order=0)

LINK_STYLE = {"color": "#1f4e8c"}
# "0 auto 1rem" = no top margin, centered, 1rem gap below each paragraph
TEXT_STYLE = {"textAlign": "center", "maxWidth": "800px", "margin": "0 auto 1rem"}
ALPHAGENOME_URL = "https://deepmind.google/science/alphagenome/"
CADD_URL = "https://cadd.gs.washington.edu/"
EVO2_URL = "https://arcinstitute.org/tools/evo"

def figure_block(filename, text):
    return html.Div(
        [
            html.Img(
                src=dash.get_asset_url(filename),
                className="img-fluid rounded d-block mx-auto mb-2",
                style={"maxWidth": "600px"},
            ),
            html.P(
                    text,
                    style={"textAlign": "center", "maxWidth": "800px", "margin": "0 auto"}
                    ),
        ],
        className="mb-4",
    )

layout = html.Div([
    html.H2("Home"),
    figure_block(
        "variants.png",
        "Variants are changes in the DNA sequence. These can be as simple as "
        "single nucleotide polymorphisms, which is shown in this image, or they "
        "can be more complex, such as insertions or deletions of larger "
        "segments of DNA.",
    ),
    figure_block(
        "mutation_effect.png",
        "Variants can have different effects on the function of genes and "
        "proteins. Some variants may have no effect, while others can lead to "
        "diseases or other health conditions.",
    ),
    figure_block(
        "prediction.png",
        "Predicting the effects of variants is a complex task that starts with "
        "analyzing the variant itself, to predicting the molecular and cellular "
        "phenotype, then using that to predict the clinical phenotypic "
        "consequences. Computational tools and algorithms are used to make "
        "these predictions.",
    ),
    html.P(
        [
            "We tested three different models: ",
            html.A("CADD", href=CADD_URL, target="_blank", style=LINK_STYLE),
            ", ",
            html.A("Evo2", href=EVO2_URL, target="_blank", style=LINK_STYLE),
            ", and ",
            html.A("AlphaGenome", href=ALPHAGENOME_URL, target="_blank", style=LINK_STYLE),
            ", against a base dataset of variants taken from MPRA with known "
            "effects.",
        ],
        style=TEXT_STYLE,
    ),
    html.P(
        "CADD is a machine learning model that combines many genomic "
        "annotations, such as conservation, regulatory, and protein-level "
        "features, into a single score of how deleterious a variant is "
        "likely to be.",
        style=TEXT_STYLE,
    ),
    html.P(
        "Evo2 is a deep learning DNA language model trained on genomes from "
        "across the tree of life that learns directly from raw sequence, and "
        "scores a variant by how much it changes the model's likelihood of "
        "the surrounding sequence.",
        style=TEXT_STYLE,
    ),
    html.P(
        "AlphaGenome is a deep learning model that also reads DNA sequence "
        "alone, but is trained to predict functional genomic signals such as "
        "gene expression and chromatin accessibility, and scores a variant "
        "by comparing predictions for the reference and alternate sequences.",
        style=TEXT_STYLE,
    ),
    html.P(
        "We chose these models because each one uses a more "
        "advanced level of machine learning than the last, and we wanted to "
        "see whether prediction performance improved as the machine "
        "learning became more sophisticated.",
        style=TEXT_STYLE,
    ),
        html.P(
        [
            "Information about how to run and test these models with your own "
            "datasets can be found in the ",
            dcc.Link("Readme", href="/readme", style=LINK_STYLE),
            " tab. The Readme tab also contains information about the models "
            "themselves, including how they were trained and what data they "
            "were trained on. This information can help you understand the "
            "strengths and limitations of each model, and how to best use them "
            "for your own research.",
        ],
        style=TEXT_STYLE,
    ),
    html.P(
        [
            "Our main hypothesis was that conservation level would have an effect on "
            "the accuracy of the predictions. Using conservation levels from the "
            "UCSC Genome Browser, we split the dataset into these groups from "
            "most to least conserved: mammals, primates, apes, and great apes. "
            "We then looked at the accuracy of each model in predicting the "
            "effects of variants in each conservation group. The results are "
            "shown in the ",
            dcc.Link("Model Analytics", href="/analytics", style=LINK_STYLE),
            " page.",
        ],
        style=TEXT_STYLE,
    ),
    html.P(
        "The results show that the accuracy of the models varies depending on "
        "the conservation level of the variants. In general, the models "
        "performed better on variants that were more conserved across species. "
        "This suggests that conservation level is an important factor to "
        "consider when predicting the effects of variants.",
        style=TEXT_STYLE,
    )
])

