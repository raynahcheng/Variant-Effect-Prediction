# Variant Model Dashboard (skeleton)

Multi-page Plotly Dash app comparing CADD, AlphaGenome and Evo2.

    pip install -r requirements.txt
    python app.py   # http://localhost:8050

| Page | Path | Contents |
|---|---|---|
| Home | `/` | Clade dropdown, model selector, ROC + per-species bar plot |
| Model Analytics | `/analytics` | Clade / model / chart dropdowns |
| Readme | `/readme` | One tab per model, loads `readmes/<model>.md` |
| Dataset | `/dataset` | Filterable table, CSV download |

Mock data lives in `data.py`; swap `VARIANTS` for real data.
