"""
Independent FP64 ground-truth builders for sdf-gems verification.

Every builder is deliberately *independent* of the formula being tested: it only uses a
definition of the shape (a polyline, a vertex set, or an inside/outside predicate).

Provides:
  - polyline_distance(P, Q, closed)        : exact distance from points Q to a dense polyline P (2D or 3D)
  - polygon_sdf(P, Q)                      : signed distance to a closed 2D polygon (dense outline of any 2D shape)
  - convex_hull_sdf(V, Q)                  : exact signed distance to the convex hull of points V (polyhedra,
                                             or convex hulls of dense curve samples, e.g. an oloid from 2 circles)
  - point_triangle_distance(Q, A, B, C)    : vectorized exact point-to-triangle distance (Ericson)
  - implicit_boundary_cloud(inside, ...)   : dense boundary point cloud of a star-shaped solid by ray-cast
                                             bisection on an inside predicate
  - cloud_sdf(cloud, Q, inside)            : signed distance via KD-tree to a dense boundary cloud, with local
                                             tangent-plane refinement for near-surface points (add exact samples
                                             of sharp edges/vertices to the cloud)

Run via UV:
  uv run tools/ground_truth.py   # self-test against analytic sphere / box distances
"""

import numpy as np
from scipy.spatial import ConvexHull, cKDTree


def polyline_distance(P, Q, closed=False, k=4):
    """Exact Euclidean distance from each row of Q to the polyline through the rows of P."""
    P = np.asarray(P, dtype=np.float64)
    Q = np.asarray(Q, dtype=np.float64)
    n = len(P)
    _, idx = cKDTree(P).query(Q, k=k)
    best = np.full(len(Q), np.inf)
    for j in range(k):
        for off in (-1, 0):
            i = idx[:, j] + off
            if closed:
                a, b = P[i % n], P[(i + 1) % n]
            else:
                i = np.clip(i, 0, n - 2)
                a, b = P[i], P[i + 1]
            ab = b - a
            t = np.clip(np.sum((Q - a) * ab, 1) / np.maximum(np.sum(ab * ab, 1), 1e-300), 0.0, 1.0)
            best = np.minimum(best, np.linalg.norm(Q - (a + t[:, None] * ab), axis=1))
    return best


def polygon_sdf(P, Q):
    """Signed distance to the closed 2D polygon P (negative inside)."""
    from matplotlib.path import Path

    d = polyline_distance(P, Q, closed=True)
    return np.where(Path(P).contains_points(Q), -d, d)


def point_triangle_distance(Q, A, B, C):
    """Vectorized exact distance from points Q (N,3) to triangle ABC (Ericson, Real-Time Collision Detection)."""
    ab, ac = B - A, C - A
    ap, bp, cp = Q - A, Q - B, Q - C
    d1, d2, d3, d4, d5, d6 = ap @ ab, ap @ ac, bp @ ab, bp @ ac, cp @ ab, cp @ ac
    va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
    den = 1.0 / np.where(np.abs(va + vb + vc) > 1e-300, va + vb + vc, 1e-300)
    res = A + np.outer(vb * den, ab) + np.outer(vc * den, ac)
    # Assign Voronoi regions in REVERSE of Ericson's early-return order so higher-priority regions win.
    m = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0)
    t = (d4 - d3) / np.maximum((d4 - d3) + (d5 - d6), 1e-300)
    res[m] = (B + np.outer(t, C - B))[m]
    m = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    t = d2 / np.maximum(d2 - d6, 1e-300)
    res[m] = (A + np.outer(t, ac))[m]
    res[(d6 >= 0) & (d5 <= d6)] = C
    m = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    t = d1 / np.maximum(d1 - d3, 1e-300)
    res[m] = (A + np.outer(t, ab))[m]
    res[(d3 >= 0) & (d4 <= d3)] = B
    res[(d1 <= 0) & (d2 <= 0)] = A
    return np.linalg.norm(Q - res, axis=1)


def convex_hull_sdf(V, Q):
    """Exact signed distance to the convex hull of V: max of facet planes inside, min triangle distance outside."""
    V = np.asarray(V, dtype=np.float64)
    Q = np.asarray(Q, dtype=np.float64)
    h = ConvexHull(V)
    inner = (Q @ h.equations[:, :3].T + h.equations[:, 3]).max(1)
    outer = np.full(len(Q), np.inf)
    ext = inner > 0
    if np.any(ext):
        Qe = Q[ext]
        de = np.full(len(Qe), np.inf)
        for s in h.simplices:
            de = np.minimum(de, point_triangle_distance(Qe, *V[s]))
        outer[ext] = de
    return np.where(ext, outer, inner)


def fibonacci_directions(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5**0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.cos(phi), np.sin(th) * np.sin(phi)], 1)


def implicit_boundary_cloud(inside, n_dirs=2_000_000, center=(0, 0, 0), r_max=2.0, iters=48, chunk=250_000):
    """Boundary samples of a solid that is star-shaped w.r.t. `center`, via bisection along rays.
    `inside(X)` -> bool array. Returns (n_dirs, 3) points on the boundary (to ~r_max * 2^-iters)."""
    c = np.asarray(center, dtype=np.float64)
    U = fibonacci_directions(n_dirs)
    out = np.empty_like(U)
    for k in range(0, n_dirs, chunk):
        u = U[k : k + chunk]
        lo = np.zeros(len(u))
        hi = np.full(len(u), r_max)
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            ins = inside(c + mid[:, None] * u)
            lo = np.where(ins, mid, lo)
            hi = np.where(ins, hi, mid)
        out[k : k + chunk] = c + (0.5 * (lo + hi))[:, None] * u
    return out


def cloud_sdf(cloud, Q, inside, refine=True, k=12):
    """Signed distance from Q to a dense boundary point cloud (sign from the `inside` predicate).

    Raw nearest-sample distance errs by up to ~spacing/2 for points close to the surface.
    With `refine=True`, points within ~4 sample spacings use the distance to a PCA tangent plane
    fitted to their k nearest samples (error ~ spacing^2 * curvature on smooth patches).
    Sharp edges/vertices: add exact dense samples of them to the cloud, and treat near-edge results
    with care (a tangent plane is ill-defined there; the code falls back to the raw distance when the
    neighbourhood is not planar)."""
    Q = np.asarray(Q, dtype=np.float64)
    tree = cKDTree(cloud)
    d, idx = tree.query(Q, k=k)
    dn = d[:, 0]
    if refine:
        spacing = np.median(tree.query(cloud[:: max(1, len(cloud) // 20000)], k=2)[0][:, 1])
        near = dn < 4.0 * spacing
        if np.any(near):
            nb = cloud[idx[near]]                       # (M, k, 3)
            c = nb.mean(1)
            X = nb - c[:, None, :]
            _, sv, vt = np.linalg.svd(X, full_matrices=False)
            nrm = vt[:, -1, :]
            planar = sv[:, -1] < 0.05 * sv[:, 0]        # neighbourhood is a smooth patch
            dp = np.abs(np.sum((Q[near] - c) * nrm, 1))
            dn = dn.copy()
            dn[near] = np.where(planar, np.minimum(dp, dn[near]), dn[near])
    return np.where(inside(Q), -dn, dn)


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    Q = rng.uniform(-2, 2, (20000, 3))
    # 1. convex hull of a cube vs exact box SDF
    V = np.array([[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)], dtype=float) * 0.5
    qa = np.abs(Q) - 0.5
    box = np.linalg.norm(np.maximum(qa, 0), axis=1) + np.minimum(qa.max(1), 0)
    print("convex_hull_sdf (cube)    max err:", np.abs(convex_hull_sdf(V, Q) - box).max())
    # 2. polygon of a dense circle vs exact disk SDF
    t = np.linspace(0, 2 * np.pi, 200000, endpoint=False)
    C = np.stack([np.cos(t), np.sin(t)], 1)
    print("polygon_sdf (disk)        max err:", np.abs(polygon_sdf(C, Q[:, :2]) - (np.linalg.norm(Q[:, :2], axis=1) - 1)).max())
    # 3. implicit cloud of a sphere
    ins = lambda X: np.linalg.norm(X, axis=1) <= 1.0
    cl = implicit_boundary_cloud(ins, n_dirs=400000)
    ex = np.linalg.norm(Q, axis=1) - 1
    print("cloud_sdf (sphere, 400k)  raw     max err:", np.abs(cloud_sdf(cl, Q, ins, refine=False) - ex).max())
    print("cloud_sdf (sphere, 400k)  refined max err:", np.abs(cloud_sdf(cl, Q, ins) - ex).max())
