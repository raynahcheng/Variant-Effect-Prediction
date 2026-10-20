import dash
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, callback, dcc, html

from data import (ALL_MODELS, CLADE_COLORS, CLADES, MODEL_OPTIONS, MODELS, VARIANTS,
                  auc, filter_by_clade, kendall_tau, roc_curve)

dash.register_page(__name__, path="/analytics", name="Model Analytics", order=1)

CHARTS = ["ROC curves", "Kendall's tau vs MPRA", "Mean score by species group",
          "Model correlation"]

# Color = species group (same as the home page); line style = model, fixed per
# model so it doesn't change as the selection changes.
MODEL_DASHES = dict(zip(ALL_MODELS, ["solid", "dash", "dashdot", "dot"]))

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
                            id="an-models", options=MODEL_OPTIONS, value=MODELS, multi=True
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
    models = [m for m in (models or []) if m in MODELS] or MODELS
    df = filter_by_clade(VARIANTS, clade)
    # "All" also shows each species group separately.
    clades = CLADES if clade == "All" else [clade]

    if chart == "ROC curves":
        fig = go.Figure()
        for c in clades:
            sub = filter_by_clade(VARIANTS, c)
            y = (sub["label"] == "functional").to_numpy().astype(int)
            for m in models:
                fpr, tpr = roc_curve(y, sub[m].to_numpy())
                fig.add_scatter(x=fpr, y=tpr, mode="lines", legendgroup=c,
                                line=dict(color=CLADE_COLORS[c], dash=MODEL_DASHES[m]),
                                name=f"{c} · {m} (AUC {auc(fpr, tpr):.3f})")
        fig.add_scatter(x=[0, 1], y=[0, 1], mode="lines", name="Chance",
                        line=dict(dash="dot", color="grey"))
        fig.update_layout(xaxis_title="False positive rate", yaxis_title="True positive rate")
    elif chart == "Kendall's tau vs MPRA":
        t = pd.DataFrame([(c, m, *kendall_tau(filter_by_clade(VARIANTS, c), m))
                          for c in clades for m in models],
                         columns=["group", "model", "tau", "n"])
        fig = px.bar(t, x="model", y="tau", color="group", barmode="group",
                     color_discrete_map=CLADE_COLORS, category_orders={"group": clades},
                     text=t["tau"].map("{:.3f}".format), hover_data={"n": ":,"})
        fig.add_hline(y=0, line_color="grey")
        fig.update_layout(yaxis_title="Kendall's tau vs MPRA |Value|",
                          legend_title_text="Species group")
    elif chart == "Mean score by species group":
        t = pd.DataFrame([(c, m, filter_by_clade(VARIANTS, c)[m].mean())
                          for c in clades for m in models],
                         columns=["group", "model", "mean"])
        # Models are on different scales, so each gets its own panel and y axis.
        fig = px.bar(t, x="group", y="mean", color="group", facet_col="model",
                     color_discrete_map=CLADE_COLORS,
                     category_orders={"group": clades, "model": models},
                     text=t["mean"].map("{:.3f}".format))
        fig.update_yaxes(matches=None, showticklabels=True, title_text="")
        fig.update_yaxes(title_text="Mean score", col=1)
        fig.update_xaxes(title_text="")
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_layout(legend_title_text="Species group")
    else:
        corr = df[["abs_value", *models]].corr(method="kendall").rename(
            index={"abs_value": "MPRA |Value|"}, columns={"abs_value": "MPRA |Value|"})
        fig = px.imshow(corr, text_auto=".2f", zmin=-1, zmax=1, color_continuous_scale="RdBu_r")

    fig.update_layout(title=f"{chart} ({clade}, n={len(df):,})", margin=dict(t=60, l=40, r=10, b=40))
    return fig
