import dash
import dash_bootstrap_components as dbc
import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, callback, dcc, html

from data import CLADES, MODELS, VARIANTS, auc, filter_by_clade, roc_curve

dash.register_page(__name__, path="/", name="Home", order=0)

layout = html.Div(
    [
        html.H2("Home"),
        dbc.Row(
            [
                # Left column: clade dropdown (whiteboard: All / Great Apes / Apes / Primates / Mammals)
                dbc.Col(
                    [
                        html.Label("Clade", className="fw-bold"),
                        dcc.Dropdown(
                            id="home-clade",
                            options=CLADES,
                            value="All",
                            clearable=False,
                        ),
                    ],
                    md=3,
                ),
                # Right column: model tabs + two plots
                dbc.Col(
                    [
                        html.Label("Models", className="fw-bold"),
                        dbc.RadioItems(
                            id="home-model",
                            options=[
                                {"label": f"{i + 1}. {m}", "value": m}
                                for i, m in enumerate(MODELS)
                            ],
                            value=MODELS[0],
                            inline=True,
                            className="mb-3",
                        ),
                        dbc.Row(
                            [
                                dbc.Col(dcc.Graph(id="home-curve"), md=6),
                                dbc.Col(dcc.Graph(id="home-bars"), md=6),
                            ]
                        ),
                    ],
                    md=9,
                ),
            ]
        ),
    ]
)


@callback(
    Output("home-curve", "figure"),
    Output("home-bars", "figure"),
    Input("home-clade", "value"),
    Input("home-model", "value"),
)
def update_home(clade, model):
    df = filter_by_clade(VARIANTS, clade)
    y = (df["label"] == "pathogenic").to_numpy().astype(int)

    # Left plot: ROC curve with chance diagonal
    fpr, tpr = roc_curve(y, df[model].to_numpy())
    curve = go.Figure()
    curve.add_scatter(x=fpr, y=tpr, mode="lines", name=f"{model} (AUC {auc(fpr, tpr):.2f})")
    curve.add_scatter(
        x=[0, 1], y=[0, 1], mode="lines", name="Chance", line=dict(dash="dot", color="grey")
    )
    curve.update_layout(
        title=f"{model}: ROC curve ({clade})",
        xaxis_title="False positive rate",
        yaxis_title="True positive rate",
        margin=dict(t=50, l=40, r=10, b=40),
        legend=dict(orientation="h", y=-0.25),
    )

    # Right plot: mean score per species for the selected clade
    means = df.groupby("species")[model].mean().sort_values(ascending=False)
    bars = go.Figure(go.Bar(x=means.index, y=means.values))
    bars.add_hline(y=float(np.mean(df[model])), line_dash="dash", annotation_text="mean")
    bars.update_layout(
        title=f"{model}: mean score by species",
        yaxis_title="Mean score",
        margin=dict(t=50, l=40, r=10, b=40),
    )
    return curve, bars
