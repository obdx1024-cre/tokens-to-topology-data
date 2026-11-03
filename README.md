# From Tokens to Topology: Companion Repository

Data, taxonomy and analysis code for the CAADRIA 2027 paper
*From Tokens to Topology: A Taxonomy of Spatial Failure Modes in AI-Generated Residential Floor Plans*.

The repository contains the 300 annotated layouts analysed in the paper (79 GPT-Image-1, 79 GPT-Image-2,
76 Nano Banana Pro, 66 Nano Banana 2), the 2,708 verified errors, the full three-tier taxonomy
(5 domains, 13 categories, 57 subcategories), and the scripts that reproduce every statistic, figure and table.

## Contents

```
data/
  layouts.csv        one row per layout (300)
  annotations.csv    one row per verified error (2,708)
  annotators.csv     layout -> anonymised annotator ID
  prompts.jsonl      prompt inputs and the full prompt sent to the image model, per layout
  images/<model>/    the 300 generated floor plans (PNG)
taxonomy/
  taxonomy.csv       domain -> category -> subcategory, severity tier and error counts
code/
  common.py                     shared data loader
  01_main_analysis_rq1_rq4.py   RQ1-RQ4 statistics
  02_rq3_block_tests.py         likelihood-ratio tests for the RQ3 predictor blocks
  03_rq4_network_bootstrap.py   co-occurrence network and bootstrap intervals (Figure 5 data)
  04_annotator_robustness.py    annotator checks reported in the Limitations
  05_figures_and_tables.py      Figures 3-5 and Tables 1-2
results/             outputs of scripts 01-04 (CSV, Markdown summary)
figures/             outputs of script 05 (PNG at 4,468 px width and PDF)
requirements.txt
```

## Data dictionary

**data/layouts.csv**

| Column | Description |
|---|---|
| `sample_id` | Layout ID (`<model>_<CubiCasa5K plan number>`) |
| `model` | GPT-Image-1, GPT-Image-2, Nano Banana Pro, Nano Banana 2 |
| `model_raw` | Model folder name used in `data/images/` |
| `morphology` | Prompted footprint morphology: Rectangular, L-Shape, Asymmetrical/Irregular, Composite/Sprawling |
| `bedroom_class`, `bedroom_count` | Prompted number of bedrooms (1-5) |
| `floor_area_m2` | Prompted floor area (m²) |
| `prompted_room_count`, `log_rooms` | Number of mandatory room labels in the prompt (every labelled space, including hallways and entrances) and its natural log |
| `error_count` | Number of verified errors in the layout |
| `image` | Path to the generated plan |

**data/annotations.csv**

| Column | Description |
|---|---|
| `error_id` | Error ID |
| `sample_id`, `model` | Layout the error belongs to |
| `domain`, `category`, `subcategory` | Position in the taxonomy |
| `severity` | Critical (blocks access), Functional (impairs use) or Cosmetic (degrades legibility); fixed per category |
| `annotator_note` | The annotator's short description, verbatim |
| `bbox_x`, `bbox_y`, `bbox_w`, `bbox_h` | Bounding box as fractions of image width and height (top-left origin) |

**data/annotators.csv**: `A1`-`A4` are the four architecture students who coded most layouts (270);
`A5`-`A9` and `NR` (not recorded) are lab members who coded 30 layouts during the initial pilot phase.
All annotators were blind to model identity.

**data/prompts.jsonl**: one JSON object per layout.

| Field | Description |
|---|---|
| `cubicasa_svg` | Source plan in CubiCasa5K |
| `extracted_rooms` | Room type, area, width and length extracted from the SVG |
| `morphology`, `bedroom_class`, `topology`, `total_area_m2` | Prompt parameters |
| `inventory_labels` | Mandatory room labels (BR1, K1, LR1, ...) |
| `narrative` | Room-by-room description with dimensions |
| `full_prompt` | The complete prompt sent to the image model (style prefix, morphology and area, narrative, labels and 3D/texture suffix) |
| `adjacency_graph_mermaid` | Prompted room adjacency graph (Mermaid) |

## Reproducing the results

```bash
pip install -r requirements.txt
python code/01_main_analysis_rq1_rq4.py   # takes a few minutes (10,000 permutations, 2,000 bootstraps)
python code/02_rq3_block_tests.py
python code/03_rq4_network_bootstrap.py
python code/04_annotator_robustness.py
python code/05_figures_and_tables.py
```

All random procedures use a fixed seed (2026).

| Result in the paper | Where it comes from |
|---|---|
| Error shares and layout prevalence with Wilson CIs (RQ1, Figure 3, Table 1) | `results/rq1_shares_prevalence.csv` |
| Errors per layout by model; χ² = 273.98, V = 0.184 and within-provider contrasts (RQ2) | `results/rq2_errors_per_layout_by_model.csv`, `results/rq2_model_category_tests.csv` |
| Severity and opening shifts (Figure 4) | computed in `05_figures_and_tables.py` |
| Negative Binomial GLM (Table 2) and log-room coefficient (RQ3) | `results/rq3_negbin_logrooms_covariate.csv` |
| Joint LR χ²(6) = 13.39, morphology LR χ²(3) = 5.83, Spearman ρ = 0.15 (RQ3) | `results/rq3_block_tests.csv` |
| 78 pairwise logistic regressions, smallest q = 0.067 (RQ4) | `results/rq4_pairwise_logistic_all78.csv` |
| Network degree, betweenness and bootstrap intervals (RQ4, Figure 5) | `results/figure5_centrality_bootstrap_intervals.csv`, `results/rq4_network_centrality_bootstrap.csv` |
| Annotator balance and robustness checks (Limitations) | printed by `04_annotator_robustness.py` |

`results/results_summary.md` summarises the main statistics in one file.

## Licence

- Code: MIT (see `LICENSE`).
- Data: the prompts are derived from CubiCasa5K. Check its licence before choosing the data licence
  (see `data/LICENSE.md`).

## Citation

If you use this repository, please cite the paper (citation to be added after publication).
