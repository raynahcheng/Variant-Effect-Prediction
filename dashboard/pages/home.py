import dash
from dash import html

dash.register_page(__name__, path="/", name="Home", order=0)

# Placeholder: to be filled in with background on variant effect prediction.
layout = html.Div(
    [
        html.H2("Home"),
        html.P("Background on variant effect prediction coming soon.", className="text-muted"),
    ]
)
