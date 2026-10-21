# `sdTruncatedOctahedron` — Exact 1-Sqrt 3D Truncated Octahedron / Kelvin Cell

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdTruncatedOctahedron](../images/polyhedra3d/sdTruncatedOctahedron.png) | ![sdTruncatedOctahedron slice 1](../images/polyhedra3d/sdTruncatedOctahedron_slice.png) | ![sdTruncatedOctahedron slice 2](../images/polyhedra3d/sdTruncatedOctahedron_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Mercury Demogroup (plane bound) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://mercury.sexy/hg_sdf/](https://mercury.sexy/hg_sdf/) |
| **Category** | `polyhedra3d` |
| **GLSL File** | [`glsl/polyhedra3d/sdTruncatedOctahedron.glsl`](../../glsl/polyhedra3d/sdTruncatedOctahedron.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, max error `1.78e-15`) |

---

## Mathematical Overview

Computes the exact 3D Euclidean distance to the space-filling Kelvin cell (14 faces, 36 edges, 24 vertices; Voronoi cell of the BCC lattice) in 1 `sqrt` via $O_h$ Weyl chamber folding.

---

## GLSL Implementation

```glsl
float sdTruncatedOctahedron( vec3 p, float s )
{
    p = abs(p);
    if (p.x < p.y) p.xy = p.yx;
    if (p.y < p.z) p.yz = p.zy;
    if (p.x < p.y) p.xy = p.yx;
    float d_hex = (p.x + p.y + p.z - 1.5*s) * 0.57735027;
    if ((p.x <= s && d_hex <= 0.0) || max(2.0*p.x - p.y - p.z, p.x + p.y - 2.0*p.z) <= 1.5*s)
        return max(p.x - s, d_hex);
    float t1 = clamp(0.5*(0.5*s - p.y + p.z), 0.0, p.z);
    vec3  d1 = vec3(p.x - s, max(p.y - 0.5*s + t1, 0.0), p.z - t1);
    float t2 = max(0.5*(0.5*s - p.x + p.y), 0.0);
    vec3  d2 = vec3(p.x - s + t2, p.y - 0.5*s - t2, p.z);
    return sqrt(min(dot(d1, d1), dot(d2, d2)));
}
```

---

## Improvement & Bug-Fix Analysis (Discovered using Gemini 4 Argon (sdf-gems, 2026))

![Before vs After Comparison](../images/comparisons/polyhedra_truncated_octahedron_slice_comp.png)

- **Before (Mercury `hg_sdf`)**: Plane bound (`26.12%` error near vertices).
- **After (`sdf-gems`)**: **Exact 3D Euclidean SDF** (`max_err = 1.78e-15`, 1 `sqrt`).

### Code Diff vs. Baseline (Mercury Demogroup (plane bound))
```diff
-    p = abs(p);
-    return max(max(p.x, max(p.y, p.z)) - r, (p.x + p.y + p.z - 1.5*r) * 0.57735027);
+    // Exact 1-sqrt Oh Weyl chamber projection onto square-hex and hex-hex edges
+    ...
+    return sqrt(min(dot(d1, d1), dot(d2, d2)));
```
