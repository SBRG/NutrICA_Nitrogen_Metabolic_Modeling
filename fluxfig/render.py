"""Draw the subnetwork. Vector output with per-reaction ids for Illustrator."""
import numpy as np

from .layout import (layout_nodes, place_reactions, relax_labels,
                     satellite_positions)
from .model import compress_genes

# Palette for highlighted metabolites: distinct at a glance, colourblind-safe
# enough to carry category information on its own.
HIGHLIGHT_PALETTE = ['#1f77b4', '#2ca02c', '#17becf', '#c2185b', '#7b3294',
                     '#d95f02', '#1a5fb4', '#4d9221']

THEME = {
    'font': 'Arial',
    'rxn_label_size': 8.0,
    'met_label_size': 7.0,
    'cof_label_size': 5.2,
    'value_label_size': 6.2,
    'met_face': '#ffffff',
    'met_edge': '#3f3f3f',
    'met_radius': 0.052,
    'cof_face': '#f5f5f5',
    'cof_edge': '#bdbdbd',
    'cof_radius': 0.022,
    'rxn_marker': 0.030,
    'sat_radius': 0.30,
    'edge_color': '#8a8a8a',
    'cof_edge_color': '#cccccc',
    'min_width': 0.8,
    'max_width': 4.2,
    'cof_width': 0.5,
    'cmap': 'RdBu_r',
    'text': '#1a1a1a',
    'muted': '#767676',
    'met_outline': True,
    'rxn_label_style': 'normal',
    'edge_neutral': '#8a8a8a',
    'highlight_scale': 1.0,
}

# Panel-a style: anonymous intermediates are flat grey and unlabelled, named
# metabolites are bold colour, and reactions are labelled with italic gene
# names instead of reaction ids.
PATHWAY_THEME = dict(THEME, **{
    'met_face': '#d6d6d6',
    'met_outline': False,
    'met_radius': 0.058,
    'met_label_size': 8.5,
    'rxn_label_size': 7.5,
    'rxn_label_style': 'italic',
    'value_label_size': 6.5,
    'rxn_marker': 0.0,
    'edge_neutral': '#9a9a9a',
    'min_width': 1.0,
    'max_width': 4.5,
    'highlight_scale': 1.5,
})

THEMES = {'default': THEME, 'pathway': PATHWAY_THEME}


def _norm_positions(pos, edges):
    """Rescale so the median edge is unit length -- makes sizes portable."""
    if len(edges):
        lens = [np.linalg.norm(pos[i] - pos[j]) for i, j in edges]
        lens = [l for l in lens if l > 1e-9]
        med = float(np.median(lens)) if lens else 1.0
    else:
        med = 1.0
    return pos / (med or 1.0)


def _build_graph(reactions, merge_metabolites=True):
    nodes, edges, index = [], [], {}

    def node(key):
        if key not in index:
            index[key] = len(nodes)
            nodes.append(key)
        return index[key]

    for r in reactions:
        ri = node(('R', r.id))
        for m in r.substrates:
            edges.append((node(('M', m if merge_metabolites else r.id + '|' + m)), ri))
        for m in r.products:
            edges.append((ri, node(('M', m if merge_metabolites else r.id + '|' + m))))
    return nodes, edges, index


def draw(reactions, values=None, out='figure.pdf', title=None,
         value_label='log$_2$ fold change', vlim=None, theme=None,
         show_cofactors=False, show_values=None, met_labels=None,
         figsize=None, dpi=300, colorbar=True, seed_positions=None,
         merge_metabolites=True, no_data_color=None, spread=1.9,
         halo=False, min_sep=0.70, style='pathway', rxn_label='genes',
         highlight=None, gene_override=None):
    """Render `reactions` to a vector file.

    values: {reaction_id: number}. Edges are coloured on a diverging scale and
    widened by |value|. Reactions with no value render thin and grey.
    style: 'pathway' (flat grey intermediates, italic gene labels, named
        metabolites in bold colour) or 'default' (every node labelled).
    rxn_label: 'genes', 'id' or 'name'.
    highlight: metabolites to name and colour -- {met: colour} or a list, which
        takes colours from HIGHLIGHT_PALETTE in order.
    gene_override: {reaction_id: 'ndh'} when the GPR lists more isozymes than
        the figure should say.
    spread: multiplies node spacing; raise it when labels feel cramped.
    halo: white outline behind text. Off by default because matplotlib renders
        stroked text as PATHS, which would make the labels uneditable in
        Illustrator; the collision pass usually makes it unnecessary.
    show_cofactors: draw ATP/NADH/H+ etc. as small satellites. Off by default:
        on a figure of more than a few reactions they add dozens of tiny labels
        and bury the story.
    """
    import matplotlib
    matplotlib.use('Agg')
    # Keep text as TEXT in the vector output. Matplotlib's default converts SVG
    # text to paths and PDF text to Type 3, both of which land in Illustrator
    # as uneditable outlines -- exactly what these figures are exported for.
    matplotlib.rcParams['svg.fonttype'] = 'none'
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, FancyArrowPatch
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    th = dict(THEMES.get(style, THEME), **(theme or {}))
    values = values or {}
    gene_override = gene_override or {}

    if highlight is None:
        highlight = {}
    elif not isinstance(highlight, dict):
        highlight = {m: HIGHLIGHT_PALETTE[i % len(HIGHLIGHT_PALETTE)]
                     for i, m in enumerate(highlight)}
    if met_labels is None:
        met_labels = (style != 'pathway')
    if show_values is None:
        # in pathway style the colourbar already carries the magnitude, and
        # per-edge numbers just crowd the figure
        show_values = (style != 'pathway')

    def rxn_text(r_):
        if rxn_label == 'id':
            return r_.id
        if rxn_label == 'name':
            return r_.name or r_.id
        return (gene_override.get(r_.id)
                or compress_genes(r_.gene_names) or r_.id)
    for r in reactions:
        v = values.get(r.id)
        r.value = v if v is None or np.isfinite(v) else None

    nodes, edges, index = _build_graph(reactions, merge_metabolites)
    pos = _norm_positions(layout_nodes(nodes, edges, seed_positions), edges)
    pos = pos * spread

    mkey = (lambda r, m: m if merge_metabolites else r.id + '|' + m)
    members = {r.id: np.array([pos[index[('M', mkey(r, m))]]
                               for m in r.substrates + r.products])
               for r in reactions}
    met_idx = [i for i, k in enumerate(nodes) if k[0] == 'M']
    rpos = place_reactions({r.id: pos[index[('R', r.id)]] for r in reactions},
                           members, spread=th['sat_radius'] * 0.9,
                           min_sep=min_sep * spread,
                           avoid=pos[met_idx] if met_idx else None,
                           avoid_sep=0.42 * spread)
    for r in reactions:
        pos[index[('R', r.id)]] = rpos[r.id]

    vals = [abs(v) for v in values.values() if v is not None and np.isfinite(v)]
    lim = vlim or (max(vals) if vals else 1.0)
    norm = Normalize(-lim, lim)
    cmap = plt.get_cmap(th['cmap'])
    grey = no_data_color or th['edge_neutral']

    def edge_style(v):
        if v is None:
            return grey, th['min_width']
        frac = min(abs(v) / lim, 1.0) if lim else 0.0
        return cmap(norm(v)), th['min_width'] + frac * (th['max_width'] - th['min_width'])

    fig, ax = plt.subplots(figsize=figsize or (7.5, 6.0))
    ax.set_aspect('equal')
    ax.axis('off')

    drawn_mets, cof_nodes, labels = set(), [], []
    R = th['met_radius'] * 72  # shrink is in points

    for r in reactions:
        rp = pos[index[('R', r.id)]]
        col, wid = edge_style(r.value)
        key = (lambda m: mkey(r, m))

        sp = [pos[index[('M', key(m))]] for m in r.substrates]
        pp = [pos[index[('M', key(m))]] for m in r.products]
        axis = ((np.mean(pp, axis=0) if pp else rp) -
                (np.mean(sp, axis=0) if sp else rp))

        gid = 'rxn-%s' % r.id
        for p in sp:
            ax.add_patch(FancyArrowPatch(
                p, rp, arrowstyle='-', connectionstyle='arc3,rad=0.05',
                color=col, linewidth=wid, shrinkA=R, shrinkB=1.5,
                capstyle='round', zorder=2, gid=gid))
        for p in pp:
            ax.add_patch(FancyArrowPatch(
                rp, p, arrowstyle='-|>', connectionstyle='arc3,rad=0.05',
                color=col, linewidth=wid, shrinkA=1.5, shrinkB=R,
                mutation_scale=7 + 1.6 * wid, capstyle='round',
                zorder=2, gid=gid))

        # the reaction itself gets a glyph, so edges meet something visible
        if th['rxn_marker'] > 0:
            ax.add_patch(Circle(rp, th['rxn_marker'], facecolor=col,
                                edgecolor='#2b2b2b', linewidth=0.5, zorder=4,
                                gid='rxnmarker-%s' % r.id))

        if show_cofactors:
            for side, group in ((-1, r.cof_sub), (1, r.cof_prod)):
                sats = satellite_positions(rp, axis, len(group),
                                           th['sat_radius'], side)
                for m, sat in zip(group, sats):
                    a, b = (sat, rp) if side < 0 else (rp, sat)
                    ax.add_patch(FancyArrowPatch(
                        a, b, arrowstyle='-|>' if side > 0 else '-',
                        connectionstyle='arc3,rad=%.2f' % (0.25 * side),
                        color=th['cof_edge_color'], linewidth=th['cof_width'],
                        mutation_scale=4, shrinkA=1, shrinkB=1, zorder=1,
                        gid='cof-%s-%s' % (r.id, m)))
                    cof_nodes.append((sat, m))

        for m, p in list(zip(r.substrates, sp)) + list(zip(r.products, pp)):
            k = key(m)
            if k in drawn_mets:
                continue
            drawn_mets.add(k)
            hl = highlight.get(m)
            rad = th['met_radius'] * (th['highlight_scale'] if hl else 1.0)
            ax.add_patch(Circle(
                p, rad, facecolor=hl or th['met_face'],
                edgecolor=th['met_edge'] if th['met_outline'] else 'none',
                linewidth=0.7 if th['met_outline'] else 0.0,
                zorder=3, gid='met-%s' % m))
            if hl or met_labels:
                labels.append({'anchor': p,
                               'off': np.array([0.0, -(rad + 0.085)]),
                               'text': m, 'size': th['met_label_size'],
                               'color': hl or th['text'],
                               'weight': 'bold' if hl else 'normal',
                               'style': 'normal', 'va': 'top',
                               'gid': 'metlabel-%s' % m, 'sub': None})

    for sat, m in cof_nodes:
        ax.add_patch(Circle(sat, th['cof_radius'], facecolor=th['cof_face'],
                            edgecolor=th['cof_edge'], linewidth=0.4, zorder=3))
        labels.append({'anchor': sat, 'off': np.array([0.0, -0.055]),
                       'text': m, 'size': th['cof_label_size'],
                       'color': th['muted'], 'weight': 'normal',
                       'style': 'normal', 'va': 'top', 'gid': None,
                       'sub': None})

    for r in reactions:
        labels.append({
            'anchor': pos[index[('R', r.id)]],
            'off': np.array([0.0, th['met_radius'] + 0.10]),
            'text': rxn_text(r),
            'size': th['rxn_label_size'], 'color': th['text'],
            'weight': 'semibold' if th['rxn_label_style'] != 'italic' else 'normal',
            'style': th['rxn_label_style'], 'va': 'bottom',
            'gid': 'rxnlabel-%s' % r.id,
            'sub': ('%+.2f' % r.value) if (show_values and r.value is not None) else None})

    # Text extents are measured in data units, so the view transform has to be
    # settled BEFORE measuring -- otherwise every measured size is wrong.
    node_pts = [pos] + ([np.array([s for s, _ in cof_nodes])] if cof_nodes else [])
    node_pts = np.vstack(node_pts)
    pad0 = 0.35 + 0.08 * max(np.ptp(node_pts[:, 0]), np.ptp(node_pts[:, 1]))
    ax.set_xlim(node_pts[:, 0].min() - pad0, node_pts[:, 0].max() + pad0)
    ax.set_ylim(node_pts[:, 1].min() - pad0, node_pts[:, 1].max() + pad0)

    # One collision pass over every label, using REAL rendered text extents
    # rather than a character-count guess: draw, measure, relax, reposition.
    # Every label is centre-anchored so the measured box and the drawn text
    # always agree; `off` just says which side of the node it starts on.
    anchors = np.array([l['anchor'] for l in labels], dtype=float)
    offs = np.array([l['off'] for l in labels], dtype=float)

    fx = [pe.withStroke(linewidth=2.0, foreground='white')] if halo else None
    for l, a, o in zip(labels, anchors, offs):
        l['artist'] = ax.text(a[0] + o[0], a[1] + o[1], l['text'],
                              ha='center', va='center', fontsize=l['size'],
                              color=l['color'], family=th['font'],
                              fontweight=l['weight'],
                              fontstyle=l.get('style', 'normal'),
                              zorder=6, gid=l['gid'],
                              path_effects=fx)
        l['sub_artist'] = None
        if l['sub']:
            l['sub_artist'] = ax.text(
                a[0] + o[0], a[1] + o[1], l['sub'], ha='center', va='center',
                fontsize=th['value_label_size'], color=th['muted'],
                family=th['font'], zorder=6, path_effects=fx)

    fig.canvas.draw()
    inv = ax.transData.inverted()

    sizes, drops = [], []
    for l in labels:
        bb = l['artist'].get_window_extent().transformed(inv)
        w, h = bb.width / 2, bb.height / 2
        drop = 0.0
        if l['sub_artist'] is not None:
            sb = l['sub_artist'].get_window_extent().transformed(inv)
            drop = h + sb.height / 2 + 0.010
            h += sb.height / 2 + 0.010 / 2
        sizes.append([w + 0.014, h + 0.014])
        drops.append(drop)
    sizes = np.array(sizes)

    obstacles = np.vstack([pos] + ([np.array([s for s, _ in cof_nodes])]
                                   if cof_nodes else []))
    P = relax_labels(anchors, offs, min_dist=0.14, sizes=sizes, max_shift=0.55,
                     obstacles=obstacles,
                     obstacle_radius=(th['met_radius'] * th['highlight_scale']
                                      + 0.055), iters=160)

    for l, p, drop in zip(labels, P, drops):
        if l['sub_artist'] is None:
            l['artist'].set_position((p[0], p[1]))
        else:
            l['artist'].set_position((p[0], p[1] + drop / 2))
            l['sub_artist'].set_position((p[0], p[1] + drop / 2 - drop))

    fig.canvas.draw()
    inv = ax.transData.inverted()
    boxes = [l['artist'].get_window_extent().transformed(inv) for l in labels]
    boxes += [l['sub_artist'].get_window_extent().transformed(inv)
              for l in labels if l['sub_artist'] is not None]
    x0 = min([b.x0 for b in boxes] + [node_pts[:, 0].min()])
    x1 = max([b.x1 for b in boxes] + [node_pts[:, 0].max()])
    y0 = min([b.y0 for b in boxes] + [node_pts[:, 1].min()])
    y1 = max([b.y1 for b in boxes] + [node_pts[:, 1].max()])
    m = 0.12
    ax.set_xlim(min(ax.get_xlim()[0], x0 - m), max(ax.get_xlim()[1], x1 + m))
    ax.set_ylim(min(ax.get_ylim()[0], y0 - m), max(ax.get_ylim()[1], y1 + m))

    if title:
        ax.set_title(title, fontsize=10.5, fontweight='bold', loc='left',
                     family=th['font'], color=th['text'], pad=10)

    if colorbar and vals:
        sm = ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01, aspect=22,
                          shrink=0.62)
        cb.set_label(value_label, fontsize=8.0, family=th['font'])
        cb.ax.tick_params(labelsize=7.0, length=2.5, width=0.5)
        for t in cb.ax.get_yticklabels():
            t.set_family(th['font'])
        cb.outline.set_linewidth(0.5)

    fmt = out.rsplit('.', 1)[-1].lower()
    fig.savefig(out, format=fmt if fmt in ('pdf', 'svg', 'png', 'eps') else 'pdf',
                bbox_inches='tight', dpi=dpi)
    plt.close(fig)
    return out
