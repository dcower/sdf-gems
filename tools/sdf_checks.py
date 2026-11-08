"""
Reusable numerical checks for sdf-gems SDFs (the five mandatory checks of the sdf-verification skill,
plus a ground-truth-free self-projection test).

All checks take numpy callables (e.g. from `tools.gpu_eval.gpu_fn`, so the literal GLSL is tested):
  sd(P)  -> (N,)   signed distance
  sdg(P) -> (N,4)  (d, grad.xyz)   [for 2D use (N,3) and pass P with 2 columns]

Provides:
  - probe_points(scale, dim)            : origin, axes, symmetry planes/diagonals (the usual NaN traps)
  - check_finite(sd, P)                 : NaN/Inf count
  - compare_gt(d, gt, label)            : max/mean |d - gt|, sign mismatches, max overstep (|d| > |gt|)
  - seam_pair_test(sd, P, eps)          : C0 continuity: |d(p + eps u) - d(p - eps u)| <= 2 eps + 1e-4
                                          (pass points ON your branch seams, not only random points)
  - eikonal_fd(sd, P, h)                : central-difference |grad d| statistics
  - check_sdg(sd, sdg, P, h)            : dg.x == d, |grad| == 1, cosine vs finite differences
  - self_projection(sd, sdg, P)         : x = p - d(p) grad(p) must land ON the surface (|d(x)| ~ 0).
                                          Wrong distances or gradients miss the surface. Needs no ground truth.
  - lipschitz_pairs(sd, P)              : max |d(a) - d(b)| / |a - b| (must be <= 1)
  - surface_zero(sd, S)                 : |d| on independently generated surface samples
  - exactness_certificate(sd, sdg, P, S): Lipschitz + zero-on-surface + self-projection  =>  exact SDF (no GT needed)

NOTE: an eikonal map (|grad d| - 1) cannot reveal every wrong SDF: a max()/min() of exact distances has
|grad d| = 1 almost everywhere and can still be badly wrong. Always also compare against ground truth.
"""

import numpy as np


def probe_points(scale=1.0, dim=3, n=401):
    s = np.linspace(-1.5 * scale, 1.5 * scale, n)
    z = np.zeros_like(s)
    if dim == 2:
        pts = [np.stack([s, z], 1), np.stack([z, s], 1), np.stack([s, s], 1), np.stack([s, -s], 1), [[0.0, 0.0]]]
    else:
        pts = [np.stack(c, 1) for c in ([s, z, z], [z, s, z], [z, z, s], [s, s, z], [s, z, s], [z, s, s], [s, s, s], [s, -s, s])]
        pts.append([[0.0, 0.0, 0.0]])
    return np.concatenate([np.asarray(p, dtype=np.float64) for p in pts])


def check_finite(sd, P, label=""):
    d = np.asarray(sd(P))
    bad = int(np.sum(~np.isfinite(d)))
    print(f"[{label}] NaN/Inf: {bad} / {d.size}")
    return bad


def compare_gt(d, gt, label=""):
    d = np.asarray(d, dtype=np.float64)
    gt = np.asarray(gt, dtype=np.float64)
    e = np.abs(d - gt)
    over = np.maximum(np.abs(d) - np.abs(gt), 0.0)
    sign = int(np.sum((np.sign(d) != np.sign(gt)) & (np.abs(gt) > 1e-6)))
    w = int(np.argmax(e))
    print(f"[{label}] max|d-gt| {e.max():.3e}  mean {e.mean():.3e}  max overstep {over.max():.3e}  "
          f"sign mismatches {sign}  (worst idx {w}: d={d[w]:.7g} gt={gt[w]:.7g})")
    return e.max()


def seam_pair_test(sd, P, eps=1e-5, label="", seed=0):
    P = np.asarray(P, dtype=np.float64)
    u = np.random.default_rng(seed).normal(size=P.shape)
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    j = np.abs(np.asarray(sd(P + eps * u)) - np.asarray(sd(P - eps * u)))
    bound = 2 * eps + 1e-4
    print(f"[{label}] seam C0: max|d+ - d-| {j.max():.3e} (bound {bound:.1e})  violations {int(np.sum(j > bound))}")
    return j.max()


def _fd_grad(sd, P, h):
    P = np.asarray(P, dtype=np.float64)
    g = []
    for k in range(P.shape[1]):
        e = np.zeros(P.shape[1])
        e[k] = h
        g.append((np.asarray(sd(P + e)) - np.asarray(sd(P - e))) / (2 * h))
    return np.stack(g, 1)


def eikonal_fd(sd, P, h=1e-3, label=""):
    g = np.linalg.norm(_fd_grad(sd, P, h), axis=1)
    e = np.abs(g - 1)
    print(f"[{label}] eikonal |FD grad|-1: median {np.median(e):.2e}  99% {np.percentile(e, 99):.2e}  "
          f"frac within 1e-2: {(e < 1e-2).mean():.4f}")
    return e


def check_sdg(sd, sdg, P, h=1e-3, label=""):
    dg = np.asarray(sdg(P), dtype=np.float64)
    d = np.asarray(sd(P), dtype=np.float64)
    dim = P.shape[1]
    g = dg[:, 1 : 1 + dim]
    fd = _fd_grad(sd, P, h)
    fdn = np.linalg.norm(fd, axis=1)
    cos = np.sum(fd * g, 1) / np.maximum(fdn * np.linalg.norm(g, axis=1), 1e-12)
    smooth = np.abs(fdn - 1) < 1e-2  # FD stencil not straddling a medial seam
    print(f"[{label}] sdg: |dg.x-d| {np.abs(dg[:, 0] - d).max():.2e}  ||g|-1| {np.abs(np.linalg.norm(g, axis=1) - 1).max():.2e}  "
          f"cos>0.999 (smooth pts) {(cos[smooth] > 0.999).mean():.4f}  NaN {int(np.sum(~np.isfinite(dg)))}")
    return cos


def self_projection(sd, sdg, P, label=""):
    P = np.asarray(P, dtype=np.float64)
    dg = np.asarray(sdg(P), dtype=np.float64)
    dim = P.shape[1]
    X = P - dg[:, :1] * dg[:, 1 : 1 + dim]
    r = np.abs(np.asarray(sd(X)))
    print(f"[{label}] self-projection: max |d(p - d grad)| {r.max():.3e}  median {np.median(r):.2e}")
    return X, r


def lipschitz_pairs(sd, P, sigma=0.05, label="", seed=1):
    """max |d(a) - d(b)| / |a - b| over pairs (a, b = a + N(0, sigma)). Must be <= 1 (+ FP32 noise)."""
    P = np.asarray(P, dtype=np.float64)
    Q = P + np.random.default_rng(seed).normal(0.0, sigma, P.shape)
    lip = np.abs(np.asarray(sd(P)) - np.asarray(sd(Q))) / np.linalg.norm(P - Q, axis=1)
    print(f"[{label}] Lipschitz: max {lip.max():.6f}  (pairs {len(P)}, sigma {sigma})")
    return lip.max()


def surface_zero(sd, S, label=""):
    """max |d| on surface samples S generated INDEPENDENTLY from the shape definition (not from sd itself)."""
    r = np.abs(np.asarray(sd(S)))
    print(f"[{label}] |d| on independent surface samples: max {r.max():.3e}")
    return r.max()


def exactness_certificate(sd, sdg, P, S, label=""):
    """Ground-truth-free exactness test. If, with S independent samples of the true surface:
      (1) d is 1-Lipschitz,
      (2) d = 0 on S (and the sign is right),
      (3) x = p - d(p) grad(p) lands on the zero set for every p,
    then |d(p)| <= dist(p, surface) by (1)+(2), and dist(p, surface) <= |p - x| = |d(p)| by (3), so d is the exact SDF
    (up to the reported tolerances). Cheap enough to run on every primitive."""
    lip = lipschitz_pairs(sd, P, label=label)
    z = surface_zero(sd, S, label=label)
    _, r = self_projection(sd, sdg, P, label=label)
    return {"lipschitz": lip, "surface_zero": z, "self_projection": float(r.max())}
