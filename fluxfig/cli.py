"""Command line entry point: python -m fluxfig ..."""
import argparse
import sys


def main(argv=None):
    p = argparse.ArgumentParser(
        prog='fluxfig',
        description='Publication figures for selected reactions of a '
                    'metabolic reconstruction, optionally coloured by fold change.')
    p.add_argument('--model', required=True, help='SBML/JSON/MAT reconstruction')
    p.add_argument('--reactions', help='comma-separated reaction ids, or a file '
                                       'with one id per line')
    p.add_argument('--values', help='CSV of values; first column reaction id')
    p.add_argument('--column', help='which column of --values to plot '
                                    '(default: the first numeric one)')
    p.add_argument('--top', type=int, help='instead of --reactions, take the N '
                                           'reactions with the largest |value|')
    p.add_argument('--min-abs', type=float, help='drop reactions with |value| below this')
    p.add_argument('--out', default='fluxfig.pdf', help='.pdf, .svg, .png or .eps')
    p.add_argument('--title')
    p.add_argument('--value-label', default='log$_2$ fold change')
    p.add_argument('--vlim', type=float, help='symmetric colour limit')
    p.add_argument('--figsize', default='7.5x6', help='inches, e.g. 9x7')
    p.add_argument('--dpi', type=int, default=300)
    p.add_argument('--cofactors', action='store_true',
                   help='draw ATP/NADH/H+ etc. as small satellite nodes')
    p.add_argument('--values-on-edges', action='store_true',
                   help='print the numeric value under each reaction label')
    p.add_argument('--style', default='pathway', choices=['pathway', 'default'])
    p.add_argument('--label', default='genes', choices=['genes', 'id', 'name'],
                   dest='rxn_label', help='what to write on each reaction')
    p.add_argument('--highlight', help='comma-separated metabolites to name and '
                                       'colour, e.g. glu__L,gln__L,nh4')
    p.add_argument('--spread', type=float, default=1.9,
                   help='node spacing; raise if labels feel cramped')
    p.add_argument('--no-met-labels', action='store_true')
    p.add_argument('--no-colorbar', action='store_true')
    p.add_argument('--split-metabolites', action='store_true',
                   help='do not merge shared metabolites into one node')
    a = p.parse_args(argv)

    from . import draw, from_cobra, load_model

    values = {}
    if a.values:
        import pandas as pd
        df = pd.read_csv(a.values, index_col=0)
        col = a.column
        if col is None:
            num = df.select_dtypes('number').columns
            if not len(num):
                p.error('no numeric column in %s' % a.values)
            col = num[0]
        elif col not in df.columns:
            p.error('column %r not in %s (have: %s)'
                    % (col, a.values, ', '.join(map(str, df.columns))))
        values = df[col].dropna().to_dict()

    if a.reactions:
        try:
            with open(a.reactions) as fh:
                ids = [l.strip() for l in fh if l.strip()]
        except OSError:
            ids = [x.strip() for x in a.reactions.split(',') if x.strip()]
    elif values:
        ids = sorted(values, key=lambda k: -abs(values[k]))
    else:
        p.error('need --reactions or --values')

    if a.min_abs is not None:
        ids = [i for i in ids if abs(values.get(i, 0.0)) >= a.min_abs]
    if a.top:
        ids = sorted(ids, key=lambda k: -abs(values.get(k, 0.0)))[:a.top]
    if not ids:
        p.error('no reactions left to draw after filtering')

    model = load_model(a.model)
    rxns, missing = from_cobra(model, ids)
    if missing:
        print('not in model, skipped: %s' % ', '.join(missing), file=sys.stderr)
    if not rxns:
        p.error('none of the requested reactions are in the model')

    w, h = (float(x) for x in a.figsize.lower().split('x'))
    out = draw(rxns, values={r.id: values.get(r.id) for r in rxns},
               out=a.out, title=a.title, value_label=a.value_label,
               vlim=a.vlim, figsize=(w, h), dpi=a.dpi,
               show_cofactors=a.cofactors,
               show_values=a.values_on_edges,
               style=a.style, rxn_label=a.rxn_label, spread=a.spread,
               highlight=[x.strip() for x in a.highlight.split(',')] if a.highlight else None,
               met_labels=not a.no_met_labels,
               colorbar=not a.no_colorbar,
               merge_metabolites=not a.split_metabolites)
    print('wrote %s (%d reactions)' % (out, len(rxns)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
