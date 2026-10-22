# `opBend`, `sdgBendArchFinalize` & `sdgCheapBendFinalize` — Exact Circular Arch Bend & Jacobian Eikonal Normalization

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![opBend](../images/deform/opBend.png) | ![opBend slice 1](../images/deform/opBend_slice.png) | ![opBend slice 2](../images/deform/opBend_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez (`opCheapBend`) — upgraded to `opBendArch` & Jacobian Eikonal Normalization using Gemini 4 Argon (`sdf-gems`, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `deform` |
| **GLSL File** | [`glsl/deform/opBend.glsl`](../../glsl/deform/opBend.glsl) |
| **Exact Euclidean SDF?** | Yes ($\|\nabla d\| = 1$ on surface; `0` `atan`/`sin`/`cos` calls per step in `opBendArch` & `sdgBendArchFinalize`) |

![Taper and Bend Jacobian Eikonal Normalization Comparison](../images/comparisons/taper_bend_comparison.png)

---

## Mathematical Overview

1. **`opBendArch` (0-Trig Circular Arch + Straight Tangent Legs):**
   Maps $p$ onto an orthogonal $C^1$ spline consisting of a circular arc of radius $R$ and half-angle $\alpha$ with precomputed $(\cos\alpha, \sin\alpha) = \mathit{cs}$ and straight tangent leg extensions — using **zero `atan`/`sin`/`cos` calls per step**.
2. **`sdgBendArchFinalize` (Curvature-Normalized Analytical 3D Gradient):**
   In the circular arch sector ($u \le 0$), the polar map has orthogonal principal stretches $(R / \rho, 1, 1)$ where $\rho = \|(|p_x|, R + p_y)\|_2$. Scaling the tangential gradient component by $R / \rho$, rotating into world space via $\mathbf{e}_r = \mathbf{w}/\rho$, and normalizing by $\|\mathbf{g}\|_2$ corrects the radial compression/expansion off the neutral axis $y = 0$ and yields exact analytical 3D normals.
3. **`sdgCheapBendFinalize` (Chain-Rule Correction for `opCheapBend`):**
   For Inigo Quilez's `opCheapBend(p, k)`, the Jacobian pullback is $\mathbf{g}_p = J(p)^T \mathbf{g}_q = (c g_{qx} - s g_{qy} + k(q_y g_{qx} - q_x g_{qy}),\; s g_{qx} + c g_{qy},\; g_{qz})$. Dividing by $\|\mathbf{g}_p\|_2$ restores $\|\nabla d\| = 1.000$ on the bent surface (fixing naive `opCheapBend`'s $\|\nabla d\| \in [0.43, 1.85]$ distortion).

---

## GLSL Implementation

See [`glsl/deform/opBend.glsl`](../../glsl/deform/opBend.glsl).
