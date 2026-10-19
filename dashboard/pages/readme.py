from pathlib import Path

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, callback, dcc, html

from data import MODELS, README_DIR

dash.register_page(__name__, path="/readme", name="Readme", order=2)

BASE = Path(__file__).resolve().parent.parent / README_DIR

layout = html.Div(
    [
        html.H2("Readme"),
        dbc.Tabs(
            id="readme-tabs",
            active_tab=MODELS[0],
            children=[dbc.Tab(label=m, tab_id=m) for m in MODELS],
            className="mb-3",
        ),
        html.Div(id="readme-body"),
    ]
)


@callback(Output("readme-body", "children"), Input("readme-tabs", "active_tab"))
def show_readme(model):
    path = BASE / f"{model.lower()}.md"
    text = path.read_text() if path.exists() else f"# {model}\n\nNo readme yet."
    return dcc.Markdown(text, link_target="_blank")
