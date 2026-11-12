import dash
import dash_bootstrap_components as dbc
from dash import Dash, dcc, html

app = Dash(
    __name__,
    use_pages=True,
    external_stylesheets=[dbc.themes.FLATLY],
    suppress_callback_exceptions=True,
    title="Variant Model Dashboard",
)

NAV_ORDER = ["Home", "Model Analytics", "Readme", "Dataset"]


def navbar():
    pages = sorted(
        dash.page_registry.values(), key=lambda p: NAV_ORDER.index(p["name"])
    )
    return dbc.NavbarSimple(
        brand="Variant Model Dashboard",
        brand_href="/",
        color="dark",
        dark=True,
        className="mb-4",
        children=[dbc.NavLink(p["name"], href=p["path"]) for p in pages],
    )


app.layout = html.Div(
    [
        dcc.Location(id="url"),
        navbar(),
        dbc.Container(dash.page_container, fluid="xl", className="pb-5"),
    ]
)

server = app.server

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8050)
