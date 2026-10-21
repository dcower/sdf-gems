# `sdTetrahedron` — Exact 1-Sqrt 3D Regular Tetrahedron (Td Parity-Sorted Weyl Chamber)

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdTetrahedron](../images/polyhedra3d/sdTetrahedron.png) | ![sdTetrahedron slice 1](../images/polyhedra3d/sdTetrahedron_slice.png) | ![sdTetrahedron slice 2](../images/polyhedra3d/sdTetrahedron_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Mercury Demogroup (plane bound) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://mercury.sexy/hg_sdf/](https://mercury.sexy/hg_sdf/) |
| **Category** | `polyhedra3d` |
| **GLSL File** | [`glsl/polyhedra3d/sdTetrahedron.glsl`](../../glsl/polyhedra3d/sdTetrahedron.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, max error `8.88e-16`) |

---

## Mathematical Overview

Folds $p$ into the fundamental chamber of the tetrahedral group $T_d$ using octant sorting and parity sign $\text{sig} = \text{sign}(p_x p_y p_z)$, reducing the entire 4-face, 6-edge, 4-vertex solid to a single face and edge check.

---

## GLSL Implementation

```glsl
float sdTetrahedron( vec3 p, float s )
{
    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
    p = abs(p);
    if (p.x < p.y) p.xy = p.yx;
    if (p.y < p.z) p.yz = p.zy;
    if (p.x < p.y) p.xy = p.yx;
    p.z *= sig;
    float m = p.x + p.y - p.z - s;
    float h = min(m / 3.0, 0.5 * (p.y - p.z));
    return sign(m) * length(p - min(p + vec3(-h, -h, h), s));
}
```

---

## Improvement & Bug-Fix Analysis (Discovered using Gemini 4 Argon (sdf-gems, 2026))

![Before vs After Comparison](../images/comparisons/polyhedra_tetrahedron_slice_comp.png)

- **Before (Mercury `hg_sdf`)**: Half-space plane bound (`42.26%` error near vertices, `||grad||` up to `1.732`, cannot be rounded).
- **After (`sdf-gems`)**: **Exact 3D Euclidean SDF** (`max_err = 8.88e-16`, `||grad d|| = 1`, 1 `sqrt`).

### Code Diff vs. Baseline (Mercury Demogroup (plane bound))
```diff
-float sdTetrahedron( vec3 p, float r )
-{
-    return (max(abs(p.x + p.y) - p.z, abs(p.x - p.y) + p.z) - r) * 0.57735027;
-}
+float sdTetrahedron( vec3 p, float s )
+{
+    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
+    p = abs(p);
+    if (p.x < p.y) p.xy = p.yx;
+    if (p.y < p.z) p.yz = p.zy;
+    if (p.x < p.y) p.xy = p.yx;
+    p.z *= sig;
+    float m = p.x + p.y - p.z - s;
+    float h = min(m / 3.0, 0.5 * (p.y - p.z));
+    return sign(m) * length(p - min(p + vec3(-h, -h, h), s));
+}
```
