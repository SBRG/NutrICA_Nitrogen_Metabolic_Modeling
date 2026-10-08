# fluxfig

Publication figures for selected reactions of a metabolic reconstruction,
optionally coloured by fold change.

Draws in the style of a hand-made pathway panel: anonymous intermediates are
flat grey circles, named metabolites are bold colour, and reactions are
labelled with italic gene names. Layout is computed from the stoichiometry, so
any reaction set can be drawn without a hand-built map.

Output is PDF / SVG / PNG / EPS with **live text** and per-reaction ids, ready
to open in Illustrator.

## Install

Nothing beyond `cobra`, `numpy`, `matplotlib`, `pandas`.

    conda activate escher_env      # already has all four

## Command line

    python -m fluxfig --model iML1515.json \
        --reactions GLNS,GLUN,GLUDy,NADH5,NADH16pp,CYTBO3_4pp,ATPS4rpp \
        --values group_log2fc_rib.csv --column D \
        --highlight glu__L,gln__L,akg,nh4,nadh,nad,q8,q8h2 \
        --title "Group D: ribosome-constrained vs optimal" \
        --out fig.pdf --figsize 8x5

Pick reactions by value instead of naming them:

    python -m fluxfig --model iML1515.json --values group_log2fc_rib.csv \
        --column C --top 12 --min-abs 1.0 --out top12.pdf

Useful flags: `--label {genes,id,name}`, `--style {pathway,default}`,
`--cofactors`, `--values-on-edges`, `--spread`, `--vlim`, `--no-colorbar`,
`--split-metabolites`.

## Python

```python
import fluxfig

model = fluxfig.load_model('iML1515.json')
rxns, missing = fluxfig.from_cobra(model, ['GLNS', 'GLUN', 'NADH5', 'NADH16pp'])

fluxfig.draw(
    rxns,
    values={'GLNS': 2.81, 'GLUN': 6.41, 'NADH5': 7.07, 'NADH16pp': -5.20},
    highlight={'glu__L': '#c2185b', 'gln__L': '#c2185b', 'nh4': '#1f77b4'},
    gene_override={'NADH5': 'ndh', 'NADH16pp': 'nuoA-N'},
    title='Group D: ribosome-constrained vs optimal',
    out='fig.pdf',
)
```

`values` is any per-reaction number — log2 fold change, flux difference,
anything. Edges take colour from a diverging scale centred on zero and width
from |value|; reactions with no value stay thin and grey.

## How it decides what to draw

* **Cofactors** (`DEFAULT_COFACTORS`) become small satellites instead of shared
  nodes, so ATP and H+ do not tie every reaction to every other one.
* A cofactor that **bridges two or more** of the selected reactions is promoted
  to a real node, because that is what links them — `q8`/`q8h2` across the
  respiratory chain, `nh4` across the GS/GOGAT reactions. The ubiquitous ones
  (`h`, `h2o`, `pi`, `atp`, `adp`, ...) never are; see `NEVER_PROMOTE`.
* A reaction written purely in cofactors (`NADH5`: nadh + h + q8 -> nad + q8h2)
  still gets a backbone, so nothing is ever drawn floating.

## Layout

Classical MDS on graph distances, refined by stress majorization, in numpy
only. Deterministic: the same input always gives the same figure. Components
are scaled to equal crowding and packed in a grid, reaction markers are snapped
onto their metabolites and fanned apart when parallel (`NADH5` vs `NADH16pp`),
and labels are placed by measuring real rendered text extents and relaxing
collisions against each other and against the nodes.

Seed known positions with `seed_positions={('M', 'glu__L'): (0, 0)}` to pin
part of a layout.

## Illustrator notes

* Text stays text: `svg.fonttype='none'`, `pdf.fonttype=42`, Arial throughout.
* Every reaction's edges are one group, `id="rxn-GLNS"`; metabolites are
  `id="met-glu__L"`, labels `id="rxnlabel-..."` / `id="metlabel-..."`.
* `halo=True` adds a white outline behind text, but matplotlib renders stroked
  text as paths, which makes those labels uneditable. Off by default.
* Label-vs-edge overlaps are not modelled (only label-vs-label and
  label-vs-node), so an occasional label sits on a line; nudge it in Illustrator.
