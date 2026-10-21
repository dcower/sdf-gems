# `sdDodecahedron` — Exact 1-Sqrt 3D Regular Dodecahedron (Ih Chamber + Axis-Aligned Edge)

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdDodecahedron](../images/polyhedra3d/sdDodecahedron.png) | ![sdDodecahedron slice 1](../images/polyhedra3d/sdDodecahedron_slice.png) | ![sdDodecahedron slice 2](../images/polyhedra3d/sdDodecahedron_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Mercury Demogroup (plane bound) — upgraded/fixed by Dillon Cower (sdf-gems, 2026)** |
| **Source** | [https://mercury.sexy/hg_sdf/](https://mercury.sexy/hg_sdf/) |
| **Category** | `polyhedra3d` |
| **GLSL File** | [`glsl/polyhedra3d/sdDodecahedron.glsl`](../../glsl/polyhedra3d/sdDodecahedron.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, max error `2.22e-15`) |

---

## Mathematical Overview

After 3 reflections into the icosahedral $I_h$ chamber, the active pentagon edge is parallel to the $y$-axis at $(x, z) = (0, 1.17557 r)$, reducing the exact 12-face, 30-edge, 20-vertex Euclidean Dodecahedron SDF to a **1-line axis-aligned segment check**!

---

## GLSL Implementation

```glsl
float sdDodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for (int i = 0; i < 3; i++) {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);
    const vec2 n = vec2(0.52573111, 0.85065081);
    float d = dot(p.xz, n) - r;
    return (p.x >= d * n.x) ? d : length(vec3(p.x, max(p.y - 0.44902798 * r, 0.0), p.z - 1.17557050 * r));
}
```

---

## Improvement & Bug-Fix Analysis (Dillon Cower (sdf-gems, 2026))

![Before vs After Comparison](../images/comparisons/polyhedra_dodecahedron_slice_comp.png)

- **Before (Mercury `hg_sdf`)**: Plane bound (`20.53%` error near vertices).
- **After (`sdf-gems`)**: **Exact 3D Euclidean SDF** (`max_err = 2.22e-15`, 1 `sqrt`).

### Code Diff vs. Baseline (Mercury Demogroup (plane bound))
```diff
-    p = abs(p);
-    const vec2 n = vec2(0.52573111, 0.85065081);
-    return max(dot(p.xy, n), max(dot(p.yz, n), dot(p.zx, n))) - r;
+    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
+    for (int i = 0; i < 3; i++) {
+        p = abs(p);
+        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
+    }
+    p.xy = abs(p.xy);
+    const vec2 n = vec2(0.52573111, 0.85065081);
+    float d = dot(p.xz, n) - r;
+    return (p.x >= d * n.x) ? d : length(vec3(p.x, max(p.y - 0.44902798 * r, 0.0), p.z - 1.17557050 * r));
```
