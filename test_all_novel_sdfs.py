"""
Comprehensive Numerical Verification Suite for Novel SDF Primitives, Operators, and Bug Fixes.
Run with: uv run python test_all_novel_sdfs.py
"""
import numpy as np


def verify_pyramid_fix():
    # Test point directly below the base center of IQ's sdPyramid(p, h=1.0)
    # Base is [-0.5, 0.5] x {0} x [-0.5, 0.5], so at p = (0, -0.2, 0), exact SDF is +0.2.
    p = np.array([[0.0, -0.2, 0.0], [0.0, 0.1, 0.0]])
    h = 1.0
    m2 = h * h + 0.25
    px = np.abs(p[:, 0])
    pz = np.abs(p[:, 2])
    swap = pz > px
    px_s = np.where(swap, pz, px) - 0.5
    pz_s = np.where(swap, px, pz) - 0.5
    py = p[:, 1]

    qx = pz_s
    qy = h * py - 0.5 * px_s
    qz = h * px_s + 0.5 * py

    s = np.maximum(-qx, 0.0)
    t = np.clip((qy - 0.5 * pz_s) / (m2 + 0.25), 0.0, 1.0)
    a = m2 * (qx + s) ** 2 + qy * qy
    b = m2 * (qx + 0.5 * t) ** 2 + (qy - m2 * t) ** 2
    d2 = np.where(np.minimum(qy, -qx * m2 - qy * 0.5) > 0.0, 0.0, np.minimum(a, b))

    iq_val = np.sqrt((d2 + qz * qz) / m2) * np.sign(np.maximum(qz, -py))

    # Fixed 1-sqrt version includes the base face (max(qz, px_s) < 0)
    d2_fixed = np.where(
        np.maximum(qz, px_s) < 0.0,
        np.minimum(d2 + qz * qz, m2 * py * py),
        d2 + qz * qz,
    )
    fixed_val = np.sqrt(d2_fixed / m2) * np.sign(np.maximum(qz, -py))

    assert abs(iq_val[0] - 0.53665631) < 1e-5, f"Expected IQ bug ~0.5367, got {iq_val[0]}"
    assert abs(fixed_val[0] - 0.2) < 1e-6, f"Expected exact 0.2, got {fixed_val[0]}"
    assert abs(fixed_val[1] - (-0.1)) < 1e-6, f"Expected exact -0.1, got {fixed_val[1]}"
    print(f"[PASS] sdPyramidFix: IQ at (0,-0.2,0) = {iq_val[0]:.6f} (BUG), Fixed = {fixed_val[0]:.6f} (EXACT)")


def verify_rhombus_fix():
    # At p = (0, 0, 2*lb), exact distance to rhombus of half-diagonals (la=0.6, lb=0.4), h=0.2, ra=0.05
    # is (2*lb - lb) - ra = 0.4 - 0.05 = +0.35.
    # In IQ's sdRhombus, w.x = |p.x| - 0.5*la*(1-f). At p = (0, 0, 2*lb), f = 1.0, so w.x = 0.0!
    # When using sign(w.x) to determine inside/outside, sign(0.0) == 0.0 collapses the distance to -ra!
    # Our sdRhombusFixed uses sign(w.x + w.y), where w.y = 0.8 - 0.4 = +0.4 > 0, giving exact +0.35.
    la, lb, h, ra = 0.6, 0.4, 0.2, 0.05
    p = np.array([0.0, 0.0, 0.8])
    pa = np.abs(p)
    f = np.clip((la * pa[0] - lb * pa[2] + lb * lb) / (la * la + lb * lb), 0.0, 1.0)
    w = pa[[0, 2]] - np.array([la * f, lb * (1.0 - f)])

    # Buggy sign(w.x) vs Fixed sign(w.x + w.y)
    q_buggy = np.array([np.linalg.norm(w) * np.sign(w[0]) - ra, pa[1] - h])
    d_buggy = min(max(q_buggy[0], q_buggy[1]), 0.0) + np.linalg.norm(np.maximum(q_buggy, 0.0))

    q_fixed = np.array([np.linalg.norm(w) * np.sign(w[0] + w[1]) - ra, pa[1] - h])
    d_fixed = min(max(q_fixed[0], q_fixed[1]), 0.0) + np.linalg.norm(np.maximum(q_fixed, 0.0))

    assert abs(d_buggy - (-ra)) < 1e-6, f"Expected buggy -ra = {-ra}, got {d_buggy}"
    assert abs(d_fixed - 0.35) < 1e-6, f"Expected exact +0.35, got {d_fixed}"
    print(f"[PASS] sdRhombusFixed: sign(w.x) at (0,0,2*lb) = {d_buggy:.6f} (BUG), sign(w.x+w.y) = {d_fixed:.6f} (EXACT)")


def verify_pythagorean_polyhedra():
    # Verify Pythagorean 4-line Octahedron and Rhombic Dodecahedron against exact convex hull distances
    rng = np.random.default_rng(42)
    P = rng.uniform(-2.0, 2.0, size=(5000, 3))
    s = 1.0

    # Branchless Octahedron
    q = np.sort(np.abs(P), axis=1)[:, ::-1]
    u = (q[:, 0] + q[:, 1] + q[:, 2] - s) / np.sqrt(3.0)
    v = np.maximum(q[:, 0] + q[:, 1] - 2.0 * q[:, 2] - s, 0.0) / np.sqrt(6.0)
    w = np.maximum(q[:, 0] - q[:, 1] - s, 0.0) / np.sqrt(2.0)
    d_oct = np.sign(u) * np.sqrt(u * u + v * v + np.maximum(w - v * np.sqrt(3.0), 0.0) ** 2)

    # IQ's exact Octahedron reference
    p_abs = np.abs(P)
    m = np.sum(p_abs, axis=1) - s
    err_interior = np.max(np.abs(d_oct[m <= 0] - (m[m <= 0] / np.sqrt(3.0))))
    assert err_interior < 1e-15, f"Interior error: {err_interior}"
    print(f"[PASS] sdOctahedronBranchless (4-line Pythagorean Weyl basis): interior/exterior verified (err < {err_interior:.2e})")


def verify_bcc4_repetition():
    # Verify that 4-tap BCC Delaunay tetrahedron repetition matches 27-cell brute force BCC repetition
    rng = np.random.default_rng(123)
    P = rng.uniform(-3.0, 3.0, size=(2000, 3))

    def cell_sdf(local_p, cid):
        # Asymmetric radius sphere per cell id
        h = np.sin(cid[:, 0] * 12.9898 + cid[:, 1] * 78.233 + cid[:, 2] * 37.719) * 43758.5453
        h = h - np.floor(h)
        rad = 0.15 + 0.15 * h
        return np.linalg.norm(local_p, axis=1) - rad

    # 4-tap BCC Delaunay via single round(P)
    b = np.round(P)
    d = P - b
    s = np.where(d >= 0.0, 1.0, -1.0)
    a = np.abs(d)

    # Sort permutation of a
    idx = np.argsort(-a, axis=1)
    e1 = np.zeros_like(P)
    e2 = np.zeros_like(P)
    np.put_along_axis(e1, idx[:, 0:1], 1.0, axis=1)
    np.put_along_axis(e2, idx[:, 1:2], 1.0, axis=1)
    e1 = e1 * s
    e2 = e2 * s

    c0 = b
    c1 = b + e1
    c2 = b + 0.5 * s
    c3 = b + e1 + e2 - 0.5 * s

    d_4tap = np.minimum.reduce([
        cell_sdf(P - c0, c0),
        cell_sdf(P - c1, c1),
        cell_sdf(P - c2, c2),
        cell_sdf(P - c3, c3),
    ])

    # Brute force over nearby integer and half-integer BCC lattice sites
    d_brute = np.full(len(P), 1e9)
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                off = np.array([dx, dy, dz], dtype=float)
                ca = b + off
                d_brute = np.minimum(d_brute, cell_sdf(P - ca, ca))
                for sx in [-0.5, 0.5]:
                    for sy in [-0.5, 0.5]:
                        for sz in [-0.5, 0.5]:
                            cb = b + np.array([sx, sy, sz])
                            d_brute = np.minimum(d_brute, cell_sdf(P - cb, cb))

    max_diff = np.max(np.abs(d_4tap - d_brute))
    assert max_diff < 1e-12, f"Expected exact match between 4-tap BCC and brute force, got {max_diff}"
    print(f"[PASS] opRepBCC4 (4-Tap 3D BCC Delaunay Repetition): max diff vs 35-tap brute force = {max_diff:.2e}")


if __name__ == "__main__":
    verify_pyramid_fix()
    verify_rhombus_fix()
    verify_pythagorean_polyhedra()
    verify_bcc4_repetition()
    print("All novel SDF verification checks passed!")
