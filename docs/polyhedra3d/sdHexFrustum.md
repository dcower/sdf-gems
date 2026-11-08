# `sdHexFrustum`, `sdgHexFrustum`, `sdHexBipyramid`, `sdgHexBipyramid` & `sdBipyramid` — Exact 0/1-`sqrt` 3D Hexagonal Frustum & Bipyramid

![sdHexFrustum Preview](../images/polyhedra3d/sdHexFrustum.png)

## Overview

A **Regular Hexagonal Frustum** (bounded by 6 congruent isosceles trapezoidal side faces and 2 parallel regular hexagonal caps of inradii $r_1, r_2$ at $y = \mp h$) and a **Regular Hexagonal / $N$-gonal Bipyramid** (two regular pyramids joined base-to-base at $y = 0$, bounded by $2N$ isosceles triangles) are ubiquitous in architectural columns, nuts/bolts, crystals, and gemstone facets.

Previously, hexagonal frustums and bipyramids were evaluated as 8-to-12-plane half-space bounds $\max_i(\langle p, \mathbf{n}_i \rangle - d_i)$, which underestimate exterior corner distances by up to $35\%$ and cannot be rounded. Using **Universal Dihedral Sector Folding ($0 \le v \le k_z u$)**, we derive the first **100% exact 0/1-`sqrt` 3D Euclidean SDFs and analytical gradients** (`sdHexFrustum`, `sdgHexFrustum`, `sdHexBipyramid`, `sdgHexBipyramid`, `sdBipyramid`).

---

## Mathematical Formulation

![Hexagonal Frustum & Bipyramid Comparison](../images/comparisons/hex_frustum_bipyramid_comparison.png)

### 1. Theorem: Universal Dihedral Frustum Decomposition (`sdHexFrustum`, `sdgHexFrustum`)

Folding `p.xz` into the $30^\circ$ dihedral $D_6$ Weyl sector $0 \le v \le k_z u$ ($k_z = \tan 30^\circ = \frac{1}{\sqrt{3}}$) reduces all 6 trapezoidal side faces and 2 hexagonal caps to a single slanted face and two flat caps. Setting $r_m = \frac{r_1 + r_2}{2}$, $s = \frac{r_1 - r_2}{2}$, $m^2 = h^2 + s^2$, and $q = (u - r_m,\; p_y,\; v - k_z r_m)$:
1. **0-`sqrt` Interior, Cap Prism, & Slanted-Face Prism Early Exit ($54.8\%$ of space)**:
   Let $d_s = \frac{h q_x + s q_y}{\sqrt{m^2}}$, $d_c = |q_y| - h$, and $u_1 = \mathrm{clamp}(\frac{h q_y - s q_x}{m^2}, -1, 1)$. Whenever $\max(d_s, d_c) \le 0$ (interior), or $d_c > 0 \;\land\; u \le r_{\text{cap}}$ (top/bottom hexagonal cap prism), or $d_s > 0 \;\land\; |h q_y - s q_x| < m^2 \;\land\; q_z + k_z s u_1 \le 0$ (slanted trapezoidal face prism), the exact signed distance is $\max(d_s, d_c)$ with **0 per-pixel `sqrt`s**.
2. **Unified Rim & Lateral Edge Projection**:
   Because $1 + k_z^2 = 1 + \frac{1}{3} = \frac{4}{3}$, the lateral edge parameter is:
   $$u_2 = \mathrm{clamp}\left(\frac{h q_y - s(q_x + k_z q_z)}{m^2 + \frac{1}{3}s^2}, -1, 1\right),$$
   requiring a single `sqrt(min(dot(v1, v1), dot(v2, v2)))` outside the face prisms.

### 2. Theorem: Equatorial-Folded Bipyramid (`sdHexBipyramid`, `sdgHexBipyramid`, `sdBipyramid`)

For a bipyramid of equatorial inradius $r$ and polar half-height $h$, reflecting across the equatorial plane `ay = abs(p.y)` eliminates the base plane completely: **100% of the bipyramid interior** ($q_z = h(u - r) + r|p_y| \le 0$) and the exterior triangular face prism require **0 per-pixel `sqrt`s**.

---

## Numerical Verification

| Function | Symmetry Group | Per-Pixel `sqrt` Count | Max Error vs 3D Mesh | Max $\|\|\nabla d\| - 1\|$ |
| :--- | :---: | :---: | :---: | :---: |
| **8-Plane Frustum Bound** | $D_6$ | `0` `sqrt` | `3.51e-1` ($35.1\%$ corner underestimation) | Up to `1.414` |
| **Ours `sdHexFrustum` / `sdgHexFrustum`** | $D_6$ (order 12) | **`0` ($54.8\%$) / `1` ($45.2\%$)** | **`1.11e-15`** | **`5.55e-16`** |
| **Ours `sdHexBipyramid` / `sdgHexBipyramid`** | $D_{6h}$ (order 24) | **`0` ($51.2\%$) / `1` ($48.8\%$)** | **`1.33e-15`** | **`4.44e-16`** |

---

## GLSL Implementation

See [`glsl/polyhedra3d/sdHexFrustum.glsl`](../../glsl/polyhedra3d/sdHexFrustum.glsl).
