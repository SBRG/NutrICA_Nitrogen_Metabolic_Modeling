# NutrICA Nitrogen: Metabolic Modeling

Genome-scale modeling of *E. coli* K-12 MG1655 growing on 19 different nitrogen sources.
We use an **ME-model** (iJL1678b, which models metabolism *and* protein/RNA production) to ask:

> When the cell is forced to express the gene programs we actually measured,
> how far does it fall from the best possible growth, and where do the resources go?

## The idea in plain words

For every nitrogen source the ME-model is solved two ways:

| Solution | Name in files | What it is |
|---|---|---|
| **Constrained** | `imod_constrained_*` | The protein spent on selected **iModulons** (co-regulated gene sets) is fixed to the levels seen in RNA-seq. Closer to the "real cell". |
| **Optimal** | `uptake_constrained_*` | Nutrient uptake is fixed to the constrained solution's values, the iModulon constraints are removed, and the model grows as fast as it can. The "ideal cell". |

Comparing the two shows what the measured gene expression costs the cell.

**iModulons used**

| Group | iModulons | Roughly means |
|---|---|---|
| 1 | NtrC-1, NtrC-2, NtrC-3 | nitrogen-starvation response |
| 2 | GadX-1, GadX-2 | acid-resistance response |
| 3 | ArcA-1, ArcA-2, Fnr-1, Fnr-2, Fnr-3 | low-oxygen / respiration response |

**Nitrogen-source groups**

| Group | Nitrogen sources |
|---|---|
| A | NH4Cl, Adenine, L-Cys |
| B | Cytosine, Cytidine |
| C | Putrescine, L-Pro, Adenosine, L-Asn, GlcNAc, Glycine, L-Orn, L-Asp, D-Ala, L-Ala, L-Arg |
| D | GABA, L-Ser, Guanosine |

Main findings explored in `plots.ipynb`: how much of the proteome the iModulons take up,
the growth rate given up, the extra ATP spent (largely a glutamine "futile cycle" between
glutaminase *glsA/ybaS* (GLUN) and glutamine synthetase *glnA* (GLNS)), and changes in
low-efficiency respiration (cytochrome bd-I / bd-II, NDH-2).

## Repository contents

### Notebooks

| Notebook | What it does |
|---|---|
| `proteome_constraints.ipynb` | **Main modeling notebook.** Averages RNA-seq replicates per nitrogen source, turns iModulon expression into proteome-fraction constraints, solves the ME-model for each source and writes the `fluxes/` tables (both the constrained and optimal solutions). |
| `plots.ipynb` | **Analysis and figures.** Reads `fluxes/` and makes all figures: proteome allocation, growth rates, ATP cost, the GLNS/GLUN cycle, respiration, and per-reaction heatmaps. Commented step by step. Figures are written to `figs/`. |
| `nitrogen_quality.ipynb` | Simple flux-balance (pFBA) analysis with the iML1515 metabolic model: how "good" each nitrogen source is, and how much of each nitrogen-containing building block it supports. |
| `transcriptomic_constraints.ipynb` | Earlier exploration: how many iModulon genes are in the ME-model, and first expression-based constraints. Reads some files (`ME_fluxes_nitrogen_sources.csv`, `ME_growth_rates_nitrogen_sources.csv`) that are not in this repository. |

### Data

| Path | Contents |
|---|---|
| `tpm_qc_nitrogen42_labeled.csv` | RNA-seq expression (TPM) after quality control; columns are `<nitrogen source>__<sample>`. |
| `log_tpm_qc.csv` | log-TPM expression for the full compendium (42 MB). |
| `sample_table_TPMcolumn_to_nitrogen_source.xlsx` | Links each RNA-seq sample to its nitrogen source and measured growth rate. |
| `imodulon_genes/` | Gene lists for iModulons (from iModulonDB, PRECISE) and ribosomal genes, as CSV. |
| `imod_genes/` | Gene sets for every iModulon, one pickled Python set per iModulon (`<name>.pkl`). |
| `iML1515.json` | iML1515 *E. coli* metabolic model (M-model), used by `nitrogen_quality.ipynb` and `fluxfig`. |
| `escher_maps/` | Escher pathway maps for viewing fluxes. |
| `rib_constrained_translation_fluxes.csv` | Translation fluxes from an earlier ribosome-constrained run. |

### Model outputs (`fluxes/`)

Rows are reactions (or genes, for translation tables); columns are nitrogen sources.
Fluxes are in mmol gDW⁻¹ h⁻¹, growth rates in h⁻¹.

| File | Contents |
|---|---|
| `uptake_constrained_fluxes.csv` | All reaction fluxes, **optimal** solution |
| `uptake_constrained_growth.csv` | Growth rate, optimal solution |
| `uptake_constrained_translation_fluxes.csv` | Protein production per gene, optimal solution |
| `imod_constrained_fluxes.csv` | All reaction fluxes, **constrained** solution |
| `imod_constrained_growth.csv` | Growth rate, constrained solution |
| `imod_constrained_translation_fluxes.csv` | Protein production per gene, constrained solution |
| `imod_constrained_uptakes.csv` | Glucose and nitrogen uptake rates per source |
| `imodulon_reaction_map.csv` | Which ME-model reactions each iModulon's genes catalyse, and whether they carry flux |
| `glns_range_<source>.csv` | Lowest and highest possible GLNS flux (and matching GLUN flux) at fixed growth |

### `fluxfig/`

A small tool for drawing pathway figures of chosen reactions, coloured by any per-reaction value
(for example, fold change). See [`fluxfig/README.md`](fluxfig/README.md).

## Requirements

- Python with `pandas`, `numpy`, `scipy`, `matplotlib`, `openpyxl`
- [`cobrapy`](https://github.com/opencobra/cobrapy) for the M-model notebooks
- [`COBRAme`](https://github.com/SBRG/cobrame) and [`ECOLIme`](https://github.com/SBRG/ecolime),
  plus a solver they support (e.g. qMINOS / SoPlex), for the ME-model notebooks.
  The ME-model is loaded from `ecolime/me_models/iJL1678b.pickle` inside the installed `ecolime` package.

`plots.ipynb` loads the ME-model for the stacked proteome figure (protein weights), the ATP sections
and the heatmaps; the other figures run from the CSVs in `fluxes/`.

## How to reproduce

1. Run `proteome_constraints.ipynb` to solve the ME-model and write `fluxes/`. This is slow:
   one ME-model solve per nitrogen source, per solution type.
2. Run `plots.ipynb` to make the figures in `figs/`. Create the folder first (`mkdir figs`);
   it is not included in the repository.

## Known issues

- `plots.ipynb`, GLNS vs GLUN coupling cell: the line defining `lim` is commented out.
  Uncomment it before running that cell, or it will raise a `NameError`.
