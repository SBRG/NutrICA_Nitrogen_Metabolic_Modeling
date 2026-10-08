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

### `environment.yml`

Conda environment for the ME-model notebooks. See **Installation** below.

### `fluxfig/`

A small tool for drawing pathway figures of chosen reactions, coloured by any per-reaction value
(for example, fold change). See [`fluxfig/README.md`](fluxfig/README.md).

## Installation

### 1. Conda environment (ME-model notebooks)

`environment.yml` recreates the environment used for `proteome_constraints.ipynb`, `plots.ipynb`
and `transcriptomic_constraints.ipynb`: Python 3.6, cobra 0.5.11 and
[COBRAme](https://github.com/SBRG/cobrame), pinned to the versions the results were made with.

```bash
conda env create -f environment.yml
conda activate cobrame_env
```

### 2. ECOLIme and the *E. coli* ME-model

[ECOLIme](https://github.com/SBRG/ecolime) holds the *E. coli* data and builds the ME-model
(iJL1678b). Install it as an editable clone, because the notebooks load the built model from
inside the package folder (`ecolime/me_models/iJL1678b.pickle`).

```bash
git clone https://github.com/SBRG/ecolime.git
cd ecolime
git checkout 1d7116b32f474dd5adbf9c72e67b22654304469d   # v0.0.9, the version used here
```

One-line fix needed with current `xlrd` (which no longer reads `.xlsx`): in
`ecolime/corrections.py`, function `correct_reaction_stoichiometries`, change

```python
df = pd.read_excel(file_name, index_col=0)
```
to
```python
df = pd.read_excel(file_name, index_col=0, engine='openpyxl')
```

Then install it and build the model (this writes `iJL1678b.pickle`, and takes a while):

```bash
pip install -e .
cd ecolime
python build_me_model.py
```

### 3. SoPlex solver

The ME-model is solved with SoPlex 3.1.1 through
[soplex_cython](https://github.com/SBRG/soplex_cython). SoPlex can't be installed with conda or pip,
because its source must be downloaded from ZIB under their academic licence.

```bash
sudo apt-get install libgmp-dev                       # macOS: brew install gmp
git clone https://github.com/SBRG/soplex_cython.git
cd soplex_cython
# download soplex-3.1.1.tgz from https://soplex.zib.de and put it in this folder
pip install .
```

Check with `python -c "import soplex"`.

### Other notebooks and tools

- `nitrogen_quality.ipynb` uses the newer cobrapy API (`cobra.flux_analysis.pfba`), so it does
  **not** run in `cobrame_env`. It was run in a separate Python 3.9 environment with cobra 0.29.1:
  `conda create -n cobrapy python=3.9 && conda activate cobrapy && pip install cobra==0.29.1 pandas matplotlib scipy jupyter`.
- `fluxfig/` needs only a recent `cobra`, `numpy`, `matplotlib` and `pandas`.

### What needs the ME-model

`plots.ipynb` loads the ME-model for the stacked proteome figure (protein weights), the ATP sections
and the heatmaps; the other figures run from the CSVs in `fluxes/`. Solving the model
(`proteome_constraints.ipynb`) also needs SoPlex.

## How to reproduce

1. Set up the environment (see **Installation**) and run `proteome_constraints.ipynb` to solve the ME-model and write `fluxes/`. This is slow:
   one ME-model solve per nitrogen source, per solution type.
2. Run `plots.ipynb` to make the figures in `figs/`. Create the folder first (`mkdir figs`);
   it is not included in the repository.

## Known issues

- `plots.ipynb`, GLNS vs GLUN coupling cell: the line defining `lim` is commented out.
  Uncomment it before running that cell, or it will raise a `NameError`.
