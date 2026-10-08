"""Pull a drawable subnetwork out of a genome-scale reconstruction."""
import re

# Metabolites that connect everything to everything. Drawn as small satellites
# on their reaction instead of as shared hubs, or the figure turns to spaghetti.
DEFAULT_COFACTORS = {
    'h', 'h2o', 'atp', 'adp', 'amp', 'pi', 'ppi', 'nad', 'nadh', 'nadp',
    'nadph', 'co2', 'o2', 'nh4', 'coa', 'q8', 'q8h2', 'mqn8', 'mql8',
    '2dmmq8', '2dmmql8', 'fad', 'fadh2', 'fmn', 'fmnh2', 'gtp', 'gdp', 'utp',
    'udp', 'ctp', 'cdp', 'so4', 'h2s', 'hco3', 'na1', 'k', 'mg2', 'fe2',
    'fe3', 'cl', 'pyr_unused',
}

COMPARTMENT_RE = re.compile(r'_([a-z]\w?)$')


def strip_compartment(mid):
    return COMPARTMENT_RE.sub('', mid)


class Reaction:
    """Minimal reaction record: what the layout and renderer actually need."""

    __slots__ = ('id', 'name', 'substrates', 'products', 'cof_sub',
                 'cof_prod', 'value', 'genes', 'gene_names')

    def __init__(self, id, name='', substrates=(), products=(),
                 cof_sub=(), cof_prod=(), genes='', gene_names=()):
        self.id = id
        self.name = name
        self.substrates = list(substrates)
        self.products = list(products)
        self.cof_sub = list(cof_sub)
        self.cof_prod = list(cof_prod)
        self.genes = genes
        self.gene_names = list(gene_names)
        self.value = None

    def __repr__(self):
        return '<Reaction %s %s -> %s>' % (self.id, self.substrates, self.products)


def from_cobra(model, reaction_ids, cofactors=None, directions=None,
               keep_compartment=False):
    """Build Reaction records for `reaction_ids` from a cobrapy model.

    directions: optional {rxn_id: +1/-1}. -1 swaps substrates and products so
    the drawn arrow follows the carried flux rather than the written equation.
    """
    cofactors = DEFAULT_COFACTORS if cofactors is None else set(cofactors)
    directions = directions or {}
    out = []
    missing = []

    for rid in reaction_ids:
        try:
            r = model.reactions.get_by_id(rid)
        except KeyError:
            missing.append(rid)
            continue

        subs, prods, cs, cp = [], [], [], []
        for met, coeff in r.metabolites.items():
            label = met.id if keep_compartment else strip_compartment(met.id)
            is_cof = strip_compartment(met.id) in cofactors
            target = (cs if is_cof else subs) if coeff < 0 else (cp if is_cof else prods)
            target.append(label)

        if directions.get(rid, 1) < 0:
            subs, prods = prods, subs
            cs, cp = cp, cs

        out.append(Reaction(rid, r.name, subs, prods, cs, cp,
                            getattr(r, 'gene_reaction_rule', ''),
                            sorted(g.name or g.id for g in getattr(r, 'genes', []))))

    promote_bridges(out)
    return out, missing


# Cofactors that are never worth drawing as a shared backbone node: they touch
# essentially every reaction and would collapse the figure into a hairball.
NEVER_PROMOTE = {'h', 'h2o', 'pi', 'ppi', 'co2', 'o2', 'hco3', 'na1', 'k',
                 'mg2', 'fe2', 'fe3', 'cl', 'so4',
                 # ATP couples every module to every other one; drawing it as a
                 # shared node drags unrelated pathways together
                 'atp', 'adp', 'amp'}


def promote_bridges(reactions, min_bridge=2, never_promote=None):
    """Promote cofactors that actually connect the selected reactions.

    A redox chain written only in cofactors (NADH5: nadh + h + q8 -> nad + q8h2)
    has no primary metabolite at all, yet q8/q8h2 are exactly what links it to
    the neighbouring reactions. Any cofactor bridging >= min_bridge of the
    selected reactions becomes a real node; the ubiquitous ones never do.
    Reactions still left with no backbone promote their own cofactors, so every
    reaction is guaranteed something to draw.
    """
    never = NEVER_PROMOTE if never_promote is None else set(never_promote)

    degree = {}
    for r in reactions:
        for m in set(r.cof_sub) | set(r.cof_prod):
            degree[m] = degree.get(m, 0) + 1

    bridging = {m for m, d in degree.items()
                if d >= min_bridge and m not in never}

    for r in reactions:
        keep_s = [m for m in r.cof_sub if m in bridging]
        keep_p = [m for m in r.cof_prod if m in bridging]
        if keep_s or keep_p:
            r.substrates += keep_s
            r.products += keep_p
            r.cof_sub = [m for m in r.cof_sub if m not in bridging]
            r.cof_prod = [m for m in r.cof_prod if m not in bridging]

    # last resort: an isolated reaction still needs a backbone
    for r in reactions:
        if not r.substrates and not r.products:
            r.substrates, r.products = r.cof_sub, r.cof_prod
            r.cof_sub, r.cof_prod = [], []
    return reactions


def load_model(path):
    """Read JSON, SBML or MAT with cobrapy, chosen by extension."""
    import cobra
    low = path.lower()
    if low.endswith('.json'):
        return cobra.io.load_json_model(path)
    if low.endswith(('.xml', '.sbml', '.xml.gz')):
        return cobra.io.read_sbml_model(path)
    if low.endswith('.mat'):
        return cobra.io.load_matlab_model(path)
    if low.endswith('.yml') or low.endswith('.yaml'):
        return cobra.io.load_yaml_model(path)
    raise ValueError('unrecognised model format: %s' % path)


def compress_genes(names, max_groups=2):
    """Render a GPR the way a pathway figure does: nuoA..nuoN -> 'nuoA-N'.

    Genes sharing a prefix collapse into one token, and the tokens join with
    'or', giving labels like 'gdhA or gltBD' or 'cyoABCD' instead of a dozen
    locus tags.
    """
    if not names:
        return ''
    groups = {}
    for n in names:
        m = re.match(r'^([a-z]{2,4}?)([A-Z]\d*)$', n)
        if m:
            groups.setdefault(m.group(1), []).append(m.group(2))
        else:
            groups.setdefault(n, [])

    tokens = []
    for prefix, sufs in groups.items():
        sufs = sorted(sufs)
        if not sufs:
            tokens.append(prefix)
        elif len(sufs) == 1:
            tokens.append(prefix + sufs[0])
        elif len(sufs) <= 4:
            tokens.append(prefix + ''.join(sufs))
        else:
            tokens.append('%s%s-%s' % (prefix, sufs[0], sufs[-1]))

    tokens.sort()
    if len(tokens) > max_groups:
        tokens = tokens[:max_groups] + ['...']
    return ' or '.join(tokens)
