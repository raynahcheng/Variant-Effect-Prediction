import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, callback, dash_table, dcc, html

from data import CLADES, MODELS, VARIANTS, filter_by_clade

dash.register_page(__name__, path="/dataset", name="Dataset", order=3)

layout = html.Div(
    [
        html.H2("Dataset"),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Label("Clade", className="fw-bold"),
                        dcc.Dropdown(id="ds-clade", options=CLADES, value="All", clearable=False),
                    ],
                    md=3,
                ),
                dbc.Col(
                    [
                        html.Label("Label", className="fw-bold"),
                        dcc.Dropdown(
                            id="ds-label",
                            options=["all", "pathogenic", "benign"],
                            value="all",
                            clearable=False,
                        ),
                    ],
                    md=3,
                ),
                dbc.Col(html.Div(id="ds-summary", className="pt-4 text-muted"), md=4),
                dbc.Col(
                    dbc.Button("Download CSV", id="ds-download-btn", color="secondary", className="mt-4"),
                    md=2,
                ),
            ],
            className="mb-3",
        ),
        dash_table.DataTable(
            id="ds-table",
            page_size=15,
            sort_action="native",
            filter_action="native",
            style_table={"overflowX": "auto"},
            style_header={"fontWeight": "bold"},
        ),
        dcc.Download(id="ds-download"),
    ]
)


def _filtered(clade, label):
    df = filter_by_clade(VARIANTS, clade)
    if label != "all":
        df = df[df["label"] == label]
    return df


@callback(
    Output("ds-table", "data"),
    Output("ds-table", "columns"),
    Output("ds-summary", "children"),
    Input("ds-clade", "value"),
    Input("ds-label", "value"),
)
def update_table(clade, label):
    df = _filtered(clade, label)
    cols = [{"name": c, "id": c} for c in df.columns]
    summary = f"{len(df):,} variants · {df['species'].nunique()} species · models: {', '.join(MODELS)}"
    return df.to_dict("records"), cols, summary


@callback(
    Output("ds-download", "data"),
    Input("ds-download-btn", "n_clicks"),
    Input("ds-clade", "value"),
    Input("ds-label", "value"),
    prevent_initial_call=True,
)
def download(n, clade, label):
    if not n or dash.ctx.triggered_id != "ds-download-btn":
        return dash.no_update
    return dcc.send_data_frame(_filtered(clade, label).to_csv, "variants.csv", index=False)
