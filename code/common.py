"""Shared loader: reads the companion data files and returns them with the column names the analysis scripts use."""
import pathlib, pandas as pd
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, RESULTS, FIGURES = ROOT/'data', ROOT/'results', ROOT/'figures'
RESULTS.mkdir(exist_ok=True); FIGURES.mkdir(exist_ok=True)

def load():
    """lab: one row per error (2,708); reg: one row per layout (300)."""
    lab = pd.read_csv(DATA/'annotations.csv').rename(columns={'subcategory':'sub_category'})
    reg = pd.read_csv(DATA/'layouts.csv').rename(columns={'error_count':'flaw_count','floor_area_m2':'floor_area','bedroom_class':'complexity'})
    return lab, reg

def load_annotators():
    """sample_id -> anonymised annotator ID (A1-A4 students; A5-A9 and NR = pilot-phase lab members)."""
    return pd.read_csv(DATA/'annotators.csv', index_col=0)
