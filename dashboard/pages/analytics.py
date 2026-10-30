import dash
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, callback, dcc, html
from plotly.subplots import make_subplots

from data import (ALL_MODELS, CLADE_COLORS, CLADES, ELEMENTS, MODEL_OPTIONS, MODELS,
                  TE_VIEWS, VARIANTS, auc, filter_by_clade, filter_by_te, kendall_tau,
                  roc_curve)

dash.register_page(__name__, path="/analytics", name="Model Analytics", order=1)

CHARTS = ["ROC curves", "Kendall's tau vs MPRA", "Mean score by species group",
          "Model correlation"]

# Color = species group (same as the home page); line style = model, fixed per
# model so it doesn't change as the selection changes. The TE split gets its
# own panels rather than another color or line style.
MODEL_DASHES = dict(zip(ALL_MODELS, ["solid", "dash", "dashdot", "dot"]))

layout = html.Div(
    [
        html.H2("Model Analytics"),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Label("Species groups", className="fw-bold"),
                        dcc.Dropdown(id="an-clades", options=CLADES, value=CLADES, multi=True),
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
                    md=2,
                ),
                dbc.Col(
                    [
                        html.Label("Transposable elements", className="fw-bold"),
                        dcc.Dropdown(id="an-te", options=list(TE_VIEWS), value="N/A",
                                     clearable=False),
                    ],
                    md=2,
                ),
                dbc.Col(
                    [
                        html.Label("Gene elements", className="fw-bold"),
                        # Empty = N/A: every element pooled.
                        dcc.Dropdown(id="an-elements", options=ELEMENTS, value=[], multi=True,
                                     placeholder="N/A (all elements)"),
                    ],
                    md=3,
                ),
                dbc.Col(
                    [
                        html.Label("Chart", className="fw-bold"),
                        dcc.Dropdown(id="an-chart", options=CHARTS, value=CHARTS[0], clearable=False),
                    ],
                    md=2,
                ),
            ],
            className="mb-4",
        ),
        dcc.Graph(id="an-graph", style={"height": "560px"}),
    ]
)


@callback(
    Output("an-graph", "figure"),
    Input("an-clades", "value"),
    Input("an-models", "value"),
    Input("an-te", "value"),
    Input("an-elements", "value"),
    Input("an-chart", "value"),
)
def update_analytics(clades, models, te, elements, chart):
    models = [m for m in (models or []) if m in MODELS] or MODELS
    clades = [c for c in CLADES if c in (clades or [])] or CLADES
    panels = TE_VIEWS.get(te, [None])
    # Gene elements filter every plot; none selected means all of them.
    base = VARIANTS[VARIANTS["Element"].isin(elements)] if elements else VARIANTS
    # Groups don't overlap, so the pooled selection is just their union.
    df = base if "Baseline" in clades else pd.concat(
        [filter_by_clade(base, c) for c in clades])
    df = pd.concat([filter_by_te(df, p) for p in panels])

    def subset(c, p):
        return filter_by_te(filter_by_clade(base, c), p)

    if chart == "ROC curves":
        fig = make_subplots(rows=1, cols=len(panels), subplot_titles=panels[0] and panels,
                            shared_yaxes=True)
        for col, p in enumerate(panels, start=1):
            for c in clades:
                sub = subset(c, p)
                y = (sub["label"] == "functional").to_numpy().astype(int)
                # No ROC without both functional and neutral variants (e.g. Great
                # Apes TE has only a handful of variants).
                if len(set(y)) < 2:
                    continue
                for m in models:
                    fpr, tpr = roc_curve(y, sub[m].to_numpy())
                    fig.add_scatter(x=fpr, y=tpr, mode="lines", legendgroup=c,
                                    line=dict(color=CLADE_COLORS[c], dash=MODEL_DASHES[m]),
                                    name=f"{c} · {m}" + (f" · {p}" if p else "")
                                    + f" (AUC {auc(fpr, tpr):.3f})", row=1, col=col)
            fig.add_scatter(x=[0, 1], y=[0, 1], mode="lines", name="Chance",
                            line=dict(dash="dot", color="grey"), showlegend=col == 1,
                            row=1, col=col)
            fig.update_xaxes(title_text="False positive rate", row=1, col=col)
        fig.update_yaxes(title_text="True positive rate", row=1, col=1)
    elif chart == "Kendall's tau vs MPRA":
        t = pd.DataFrame([(p or "All", c, m, *kendall_tau(subset(c, p), m))
                          for p in panels for c in clades for m in models],
                         columns=["te", "group", "model", "tau", "n"])
        fig = px.bar(t, x="model", y="tau", color="group", barmode="group",
                     facet_col="te" if len(panels) > 1 else None,
                     color_discrete_map=CLADE_COLORS, category_orders={"group": clades},
                     text=t["tau"].map("{:.3f}".format), hover_data={"n": ":,"})
        fig.add_hline(y=0, line_color="grey")
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_yaxes(title_text="")
        fig.update_yaxes(title_text="Kendall's tau vs MPRA |Value|", col=1)
        fig.update_layout(legend_title_text="Species group")
    elif chart == "Mean score by species group":
        t = pd.DataFrame([(p or "All", c, m, subset(c, p)[m].mean())
                          for p in panels for c in clades for m in models],
                         columns=["te", "group", "model", "mean"])
        # Models are on different scales, so each gets its own panel and y axis;
        # the TE split, if any, adds a row per subset.
        fig = px.bar(t, x="group", y="mean", color="group", facet_col="model",
                     facet_row="te" if len(panels) > 1 else None,
                     color_discrete_map=CLADE_COLORS,
                     category_orders={"group": clades, "model": models, "te": panels},
                     text=t["mean"].map("{:.3f}".format))
        fig.update_yaxes(matches=None, showticklabels=True, title_text="")
        fig.update_yaxes(title_text="Mean score", col=1)
        fig.update_xaxes(title_text="")
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_layout(legend_title_text="Species group")
    else:
        fig = make_subplots(rows=1, cols=len(panels), subplot_titles=panels[0] and panels)
        for col, p in enumerate(panels, start=1):
            corr = filter_by_te(df, p)[["abs_value", *models]].corr(method="kendall").rename(
                index={"abs_value": "MPRA |Value|"}, columns={"abs_value": "MPRA |Value|"})
            fig.add_trace(go.Heatmap(z=corr.values, x=corr.columns, y=corr.index,
                                     coloraxis="coloraxis", texttemplate="%{z:.2f}"),
                          row=1, col=col)
        fig.update_yaxes(autorange="reversed")
        fig.update_layout(coloraxis=dict(cmin=-1, cmax=1, colorscale="RdBu_r"))

    te_label = "" if te == "N/A" else f", {te}"
    if not elements:
        el_label = ""
    elif len(elements) <= 3:
        el_label = f", {', '.join(elements)}"
    else:
        el_label = f", {len(elements)} elements"
    fig.update_layout(title=f"{chart} ({', '.join(clades)}{te_label}{el_label}, n={len(df):,})",
                      margin=dict(t=60, l=40, r=10, b=40))
    return fig
