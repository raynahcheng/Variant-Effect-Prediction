import dash
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, callback, dcc, html

from data import CLADES, MODELS, VARIANTS, auc, filter_by_clade, roc_curve

dash.register_page(__name__, path="/analytics", name="Model Analytics", order=1)

CHARTS = ["Score distribution", "ROC curves", "Model correlation", "Score by chromosome"]

layout = html.Div(
    [
        html.H2("Model Analytics"),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Label("Clade", className="fw-bold"),
                        dcc.Dropdown(id="an-clade", options=CLADES, value="All", clearable=False),
                    ],
                    md=3,
                ),
                dbc.Col(
                    [
                        html.Label("Models", className="fw-bold"),
                        dcc.Dropdown(
                            id="an-models", options=MODELS, value=MODELS, multi=True
                        ),
                    ],
                    md=5,
                ),
                dbc.Col(
                    [
                        html.Label("Chart", className="fw-bold"),
                        dcc.Dropdown(id="an-chart", options=CHARTS, value=CHARTS[0], clearable=False),
                    ],
                    md=4,
                ),
            ],
            className="mb-4",
        ),
        dcc.Graph(id="an-graph", style={"height": "560px"}),
    ]
)


@callback(
    Output("an-graph", "figure"),
    Input("an-clade", "value"),
    Input("an-models", "value"),
    Input("an-chart", "value"),
)
def update_analytics(clade, models, chart):
    models = models or MODELS
    df = filter_by_clade(VARIANTS, clade)

    if chart == "Score distribution":
        long = df.melt(id_vars=["label"], value_vars=models, var_name="model", value_name="score")
        fig = px.violin(long, x="model", y="score", color="label", box=True)
    elif chart == "ROC curves":
        y = (df["label"] == "pathogenic").to_numpy().astype(int)
        fig = go.Figure()
        for m in models:
            fpr, tpr = roc_curve(y, df[m].to_numpy())
            fig.add_scatter(x=fpr, y=tpr, mode="lines", name=f"{m} (AUC {auc(fpr, tpr):.2f})")
        fig.add_scatter(x=[0, 1], y=[0, 1], mode="lines", name="Chance",
                        line=dict(dash="dot", color="grey"))
        fig.update_layout(xaxis_title="False positive rate", yaxis_title="True positive rate")
    elif chart == "Model correlation":
        fig = px.imshow(df[models].corr(), text_auto=".2f", zmin=-1, zmax=1,
                        color_continuous_scale="RdBu_r")
    else:
        long = df.melt(id_vars=["chrom"], value_vars=models, var_name="model", value_name="score")
        order = [f"chr{c}" for c in range(1, 23)]
        fig = px.box(long, x="chrom", y="score", color="model",
                     category_orders={"chrom": order})

    fig.update_layout(title=f"{chart} ({clade})", margin=dict(t=60, l=40, r=10, b=40))
    return fig
