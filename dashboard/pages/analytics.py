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

# Shown above the plot for the selected chart: what it shows, then how it
# bears on the hypothesis that models predict better in more conserved elements.
# Species groups run from most conserved (Mammals) to least (Great Apes).
CHART_SUMMARIES = {
    "ROC curves": (
        "How well each model separates functional variants (significant MPRA "
        "effect, p < 1e-5) from neutral ones. Each curve sweeps the model's score "
        "threshold; the closer it bends to the top-left corner, and the higher its "
        "AUC (in the legend), the better. The dotted diagonal is chance (AUC 0.5).",
        "If conservation helps prediction, AUC should be highest for Mammals and "
        "fall off through Primates, Apes and Great Apes, with Baseline (all "
        "variants) in between. Compare the same line style (model) across colors.",
    ),
    "Kendall's tau vs MPRA": (
        "Rank agreement between each model's score and the measured MPRA effect "
        "size (|Value|). Tau of 1 means the model orders variants exactly as the "
        "MPRA does; 0 means no relationship. Hover a bar to see how many variants "
        "it is based on.",
        "Unlike the ROC, this checks whether models get the size of the effect "
        "right, not just whether there is one. If conservation helps prediction, "
        "the bars for each model should shrink from Mammals down to Great Apes.",
    ),
    "Mean score by species group": (
        "The average score each model gives to variants in each species group. "
        "Models are on different scales, so each has its own panel and y axis; "
        "compare groups within a panel, not across panels.",
        "This shows whether models see variants in more conserved elements as "
        "more damaging. Higher scores for Mammals than for Great Apes fit the "
        "hypothesis, but on their own they say nothing about accuracy, so read "
        "them alongside the ROC and Kendall's tau charts.",
    ),
    "Model correlation": (
        "Kendall's tau between every pair of models and the MPRA |Value|, pooled "
        "over the selected species groups. Red is agreement, blue is disagreement, "
        "white is none. The \"MPRA |Value|\" row shows how well each model tracks "
        "the measured effect.",
        "Select one species group at a time and compare: if conservation helps "
        "prediction, the MPRA |Value| row should be redder for Mammals than for "
        "Great Apes, and the models should also agree more with each other in "
        "conserved elements.",
    ),
}

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
        dbc.Alert(id="an-summary", color="light", className="mb-3"),
        html.Div(
            [
                dbc.Button("Download PNG", id="an-dl-png", color="secondary", size="sm",
                           className="me-2"),
                dbc.Button("Download SVG", id="an-dl-svg", color="secondary", size="sm",
                           outline=True),
            ],
            className="d-flex justify-content-end mb-2",
        ),
        # Descriptive filename for the current plot, set alongside the figure.
        dcc.Store(id="an-filename"),
        html.Div(id="an-dl-dummy", hidden=True),
        dcc.Graph(id="an-graph", style={"height": "560px"}),
    ]
)


def plot_filename(chart, clades, models, te, elements):
    """Filename describing the plot, e.g.
    ROC_curves__Mammals-Primates__AlphaGenome-CADD__TE-split__all-elements."""
    def clean(s):
        s = "".join(ch if ch.isalnum() else "-" for ch in str(s).replace("'", ""))
        return "-".join(filter(None, s.split("-")))

    def join(xs):
        return "-".join(clean(x) for x in xs)

    parts = [
        "_".join(filter(None, clean(chart).split("-"))),
        "all-groups" if set(clades) == set(CLADES) else join(clades),
        "all-models" if set(models) == set(MODELS) else join(models),
    ]
    if te != "N/A":
        te_part = clean(te)
        parts.append(te_part if te_part.upper().startswith("TE") else "TE-" + te_part)
    if not elements:
        parts.append("all-elements")
    elif len(elements) <= 3:
        parts.append(join(elements))
    else:
        parts.append(f"{len(elements)}-elements")
    return "__".join(parts)


# Export the plot exactly as it is in the browser (current zoom, hidden traces,
# size), so no server-side image engine is needed.
dash.clientside_callback(
    """
    function(nPng, nSvg, filename) {
        const trig = dash_clientside.callback_context.triggered;
        if (!trig.length || !trig[0].value) { return dash_clientside.no_update; }
        const format = trig[0].prop_id.startsWith("an-dl-svg") ? "svg" : "png";
        const gd = document.querySelector("#an-graph .js-plotly-plot");
        if (gd) {
            Plotly.downloadImage(gd, {format: format, filename: filename || "plot",
                                      scale: format === "png" ? 2 : 1});
        }
        return "";
    }
    """,
    Output("an-dl-dummy", "children"),
    Input("an-dl-png", "n_clicks"),
    Input("an-dl-svg", "n_clicks"),
    dash.State("an-filename", "data"),
    prevent_initial_call=True,
)


@callback(Output("an-summary", "children"), Input("an-chart", "value"))
def update_summary(chart):
    shows, hypothesis = CHART_SUMMARIES[chart]
    return [
        html.P([html.Strong("What this shows: "), shows], className="mb-2"),
        html.P([html.Strong("Why it matters: "), hypothesis], className="mb-0"),
    ]


@callback(
    Output("an-graph", "figure"),
    Output("an-filename", "data"),
    Output("an-graph", "config"),
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
    filename = plot_filename(chart, clades, models, te, elements)
    # The modebar camera icon uses the same name.
    config = {"toImageButtonOptions": {"filename": filename, "scale": 2}}
    return fig, filename, config
