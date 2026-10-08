"""Deterministic graph layout: classical MDS seed + stress majorization.

numpy only. Small metabolic subnetworks (5-100 reactions) lay out in
milliseconds, and the same input always gives the same figure -- unlike a
force-directed simulation.
"""
import numpy as np


def _bfs_distances(adj, n):
    """All-pairs shortest path on an unweighted graph, as a dense matrix."""
    D = np.full((n, n), np.inf)
    for s in range(n):
        D[s, s] = 0.0
        frontier, d = [s], 0
        seen = {s}
        while frontier:
            d += 1
            nxt = []
            for u in frontier:
                for v in adj[u]:
                    if v not in seen:
                        seen.add(v)
                        D[s, v] = d
                        nxt.append(v)
            frontier = nxt
    return D


def _components(adj, n):
    seen, comps = set(), []
    for s in range(n):
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        comps.append(sorted(comp))
    return comps


def _classical_mds(D):
    """Torgerson scaling -- a deterministic starting point."""
    n = len(D)
    if n == 1:
        return np.zeros((1, 2))
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J.dot(D ** 2).dot(J)
    vals, vecs = np.linalg.eigh(B)
    idx = np.argsort(vals)[::-1][:2]
    L = np.sqrt(np.maximum(vals[idx], 0.0))
    return vecs[:, idx] * L


def _stress_majorize(D, X, iters=300, tol=1e-6):
    """SMACOF: move each point to the weighted mean of its ideal positions."""
    n = len(D)
    if n < 3:
        return X
    W = np.where(D > 0, 1.0 / np.maximum(D, 1e-9) ** 2, 0.0)
    np.fill_diagonal(W, 0.0)
    Wsum = W.sum(axis=1)
    prev = None
    for _ in range(iters):
        diff = X[:, None, :] - X[None, :, :]
        dist = np.linalg.norm(diff, axis=2)
        np.fill_diagonal(dist, 1.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(dist > 1e-9, D / dist, 0.0)
        B = -W * ratio
        np.fill_diagonal(B, 0.0)
        np.fill_diagonal(B, -B.sum(axis=1))
        X = (B.dot(X)) / np.maximum(Wsum, 1e-9)[:, None]
        stress = float((W * (dist - D) ** 2).sum())
        if prev is not None and abs(prev - stress) < tol * max(prev, 1e-9):
            break
        prev = stress
    return X


def _pack(blocks, gap_frac=0.18):
    """Arrange disconnected components in a row-major grid, no overlap."""
    if len(blocks) == 1:
        return blocks
    sizes = [(float(np.ptp(b[:, 0])) or 1.0, float(np.ptp(b[:, 1])) or 1.0)
             for b in blocks]
    gap = gap_frac * max(max(w, h) for w, h in sizes)
    cols = int(np.ceil(np.sqrt(len(blocks))))
    col_w = [0.0] * cols
    row_h = []
    for i, (w, h) in enumerate(sizes):
        c, r = i % cols, i // cols
        col_w[c] = max(col_w[c], w)
        if r >= len(row_h):
            row_h.append(0.0)
        row_h[r] = max(row_h[r], h)
    xoff = np.cumsum([0.0] + [w + gap for w in col_w])
    yoff = np.cumsum([0.0] + [h + gap for h in row_h])
    placed = []
    for i, b in enumerate(blocks):
        c, r = i % cols, i // cols
        b = b - b.min(axis=0)
        b = b + np.array([xoff[c], -yoff[r]])
        placed.append(b)
    return placed


def layout_nodes(nodes, edges, seed_positions=None):
    """Position every node. `edges` are (i, j) index pairs into `nodes`."""
    n = len(nodes)
    adj = [set() for _ in range(n)]
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)

    blocks, orders = [], []
    for comp in _components(adj, n):
        idx = {g: k for k, g in enumerate(comp)}
        sub_adj = [[idx[v] for v in adj[g] if v in idx] for g in comp]
        D = _bfs_distances(sub_adj, len(comp))
        finite = D[np.isfinite(D)]
        D[~np.isfinite(D)] = (finite.max() if finite.size else 1.0) * 1.5
        X = _classical_mds(D)
        X = _stress_majorize(D, X)

        # scale so the closest pair of nodes sits at unit distance: this
        # equalises crowding across components, which matching edge length
        # does not (a 6-leaf star and a tight cycle need different scales)
        if len(X) > 1:
            d = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
            d = d[d > 1e-9]
            if d.size:
                X = X / float(np.percentile(d, 5))

        blocks.append(X)
        orders.append(comp)

    pos = np.zeros((n, 2))
    for comp, block in zip(orders, _pack(blocks)):
        pos[comp] = block

    if seed_positions:
        for k, xy in seed_positions.items():
            if k in nodes:
                pos[nodes.index(k)] = xy
    return pos


def satellite_positions(center, axis, count, radius, side=1):
    """Fan `count` cofactor nodes off one side of a reaction marker."""
    if count == 0:
        return np.zeros((0, 2))
    ax = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(ax)
    ax = ax / norm if norm > 1e-9 else np.array([1.0, 0.0])
    perp = np.array([-ax[1], ax[0]]) * side
    spread = np.linspace(-0.45, 0.45, count) if count > 1 else np.array([0.0])
    return np.array([
        center + perp * radius + ax * radius * s for s in spread
    ])


def relax_labels(anchors, offsets, min_dist, iters=80, max_shift=None,
                 sizes=None, obstacles=None, obstacle_radius=0.06):
    """Nudge labels apart so they stop overprinting (the Escher failure mode).

    sizes: optional per-label (half-width, half-height) for anisotropic
    separation -- wide text needs more horizontal room than vertical.
    max_shift: cap how far a label may drift from its anchor, so it stays
    attached to the thing it names.
    """
    P = anchors + offsets
    n = len(P)
    if n < 2:
        return P
    if sizes is None:
        sizes = np.tile(np.array([min_dist / 2, min_dist / 2]), (n, 1))
    sizes = np.asarray(sizes, dtype=float)

    for _ in range(iters):
        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                d = P[i] - P[j]
                want = sizes[i] + sizes[j]
                overlap = want - np.abs(d)
                if overlap[0] > 0 and overlap[1] > 0:
                    # push apart along whichever axis needs the smaller move
                    k = 0 if overlap[0] < overlap[1] else 1
                    sign = 1.0 if d[k] >= 0 else -1.0
                    shift = np.zeros(2)
                    shift[k] = sign * overlap[k] * 0.5
                    P[i] += shift
                    P[j] -= shift
                    moved = True
        if obstacles is not None and len(obstacles):
            for i in range(n):
                d = P[i] - obstacles
                dist = np.linalg.norm(d, axis=1)
                hit = dist < obstacle_radius
                if hit.any():
                    k = int(np.argmin(np.where(hit, dist, np.inf)))
                    v, dv = d[k], dist[k]
                    if dv < 1e-9:
                        v, dv = np.array([0.0, 1.0]), 1.0
                    P[i] += v / dv * (obstacle_radius - dv)
                    moved = True

        if max_shift is not None:
            delta = P - anchors
            dist = np.linalg.norm(delta, axis=1, keepdims=True)
            over = (dist > max_shift).ravel()
            if over.any():
                P[over] = anchors[over] + delta[over] / dist[over] * max_shift
        if not moved:
            break
    return P


def place_reactions(rxn_pos, member_positions, spread=0.26, min_sep=0.55,
                    iters=120, avoid=None, avoid_sep=None):
    """Snap each reaction marker onto its metabolites, then separate parallels.

    A reaction drawn at the centroid of the metabolites it touches makes the
    edges read as one continuous path. Reactions with an identical metabolite
    set (NADH5 vs NADH16pp) would then coincide exactly, so they are fanned out
    perpendicular to their own axis.
    """
    pos = dict(rxn_pos)
    for rid, mets in member_positions.items():
        if len(mets):
            pos[rid] = np.mean(mets, axis=0)

    groups = {}
    for rid, mets in member_positions.items():
        key = tuple(np.round(pos[rid], 6))
        groups.setdefault(key, []).append(rid)

    for key, members in groups.items():
        if len(members) < 2:
            continue
        members = sorted(members)
        base = np.array(key, dtype=float)
        axis = None
        for rid in members:
            m = member_positions[rid]
            if len(m) >= 2:
                axis = m[-1] - m[0]
                break
        if axis is None or np.linalg.norm(axis) < 1e-9:
            axis = np.array([1.0, 0.0])
        axis = axis / np.linalg.norm(axis)
        perp = np.array([-axis[1], axis[0]])
        offs = np.linspace(-1.0, 1.0, len(members))
        for rid, o in zip(members, offs):
            pos[rid] = base + perp * o * spread

    # Reactions that merely SHARE metabolites (GLNS/GLUN/GLUDy/GLUSy all sit on
    # glutamate) land on near-identical centroids and pile up. Push every pair
    # of markers apart, then pull each one back toward its own metabolites so it
    # still reads as belonging to them.
    ids = sorted(pos)
    avoid = None if avoid is None else np.asarray(avoid, dtype=float)
    asep = min_sep * 0.7 if avoid_sep is None else avoid_sep

    if len(ids) > 1 and min_sep:
        home = {rid: pos[rid].copy() for rid in ids}
        for _ in range(iters):
            moved = False
            for i in range(len(ids)):
                for j in range(i + 1, len(ids)):
                    a, b = ids[i], ids[j]
                    d = pos[a] - pos[b]
                    dist = float(np.linalg.norm(d))
                    if dist < min_sep:
                        if dist < 1e-9:
                            d, dist = np.array([1.0, 0.0]), 1.0
                        push = d / dist * (min_sep - dist) * 0.5
                        pos[a] = pos[a] + push
                        pos[b] = pos[b] - push
                        moved = True

            # a marker sitting on top of a metabolite node reads as one blob
            if avoid is not None and len(avoid):
                for rid in ids:
                    d = pos[rid] - avoid
                    dist = np.linalg.norm(d, axis=1)
                    hit = dist < asep
                    if hit.any():
                        k = int(np.argmin(np.where(hit, dist, np.inf)))
                        v, dv = d[k], dist[k]
                        if dv < 1e-9:
                            v, dv = np.array([0.0, 1.0]), 1.0
                        pos[rid] = pos[rid] + v / dv * (asep - dv)
                        moved = True

            for rid in ids:
                pos[rid] += (home[rid] - pos[rid]) * 0.10
            if not moved:
                break
    return pos
