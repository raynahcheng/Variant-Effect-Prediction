import dash
import dash_bootstrap_components as dbc
from dash import dcc, html

dash.register_page(__name__, path="/", name="Home", order=0)

LINK_STYLE = {"color": "#1f4e8c"}

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
    figure_block("variants.png", "Variants are changes in the DNA sequence. These can be as simple as single nucleotide polymorphisms, which is shown in this image, or they can be more complex, such as insertions or deletions of larger segments of DNA."),
    figure_block("mutation_effect.png", "Variants can have different effects on the function of genes and proteins. Some variants may have no effect, while others can lead to diseases or other health conditions."),
    figure_block("prediction.png", "Predicting the effects of variants is a complex task that starts with analyzing the variant itself, to predicting the molecular and cellular phenotype, then using that to predict the clinical phenotypic consequences. Computational tools and algorithms are used to make these predictions."),
    html.P(["We tested three different models: AlphaGenome, CADD, and Evo2, against a base dataset of variants taken from MPRA with known effects. We wanted to look at if conservation level had an effect on the accuracy of the predictions. Using conservation levels from the UCSC Genome Browser, we split the dataset into these groups from most to least conserved: mammals, primates, great apes, and apes. We then looked at the accuracy of each model in predicting the effects of variants in each conservation group. The results are shown in the ",
            dcc.Link("Model Analytics", href="/analytics", style=LINK_STYLE),
            " page."],
            style={"textAlign": "center", "maxWidth": "800px", "margin": "0 auto"}),
    html.P(""),
    html.P(["Information about how to run and test these models with your own datasets can be found in the ",
            dcc.Link("Readme", href="/readme", style=LINK_STYLE),
            " tab. The Readme tab also contains information about the models themselves, including how they were trained and what data they were trained on. This information can help you understand the strengths and limitations of each model, and how to best use them for your own research."],
            style={"textAlign": "center", "maxWidth": "800px", "margin": "0 auto"}),
    html.P(""),
    html.P("The results show that the accuracy of the models varies depending on the conservation level of the variants. In general, the models performed better on variants that were more conserved across species. This suggests that conservation level is an important factor to consider when predicting the effects of variants.",
            style={"textAlign": "center", "maxWidth": "800px", "margin": "0 auto"}),
])

