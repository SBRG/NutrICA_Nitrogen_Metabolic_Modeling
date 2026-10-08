"""fluxfig -- publication figures for selected reactions of a metabolic reconstruction.

    import fluxfig
    model = fluxfig.load_model('iML1515.json')
    rxns, missing = fluxfig.from_cobra(model, ['GLNS', 'GLUN', 'NADH5'])
    fluxfig.draw(rxns, values={'GLNS': 3.21, 'NADH5': 7.07}, out='fig.pdf')

Layout is computed from the stoichiometry (classical MDS + stress
majorization), so any reaction set can be drawn without a hand-made map, and
the same input always produces the same figure.
"""
from .model import (DEFAULT_COFACTORS, NEVER_PROMOTE, Reaction,
                    compress_genes, from_cobra, load_model, promote_bridges,
                    strip_compartment)
from .layout import layout_nodes
from .render import HIGHLIGHT_PALETTE, PATHWAY_THEME, THEME, THEMES, draw

__all__ = ['load_model', 'from_cobra', 'draw', 'Reaction', 'THEME',
           'DEFAULT_COFACTORS', 'NEVER_PROMOTE', 'promote_bridges',
           'layout_nodes', 'strip_compartment', 'compress_genes', 'THEMES',
           'PATHWAY_THEME', 'HIGHLIGHT_PALETTE']
__version__ = '0.1.0'
