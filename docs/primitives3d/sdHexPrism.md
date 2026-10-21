# `sdHexPrism` — Exact 3D Hexagonal Prism (1-Sqrt Fundamental Sector Optimization)

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdHexPrism](../images/primitives3d/sdHexPrism.png) | ![sdHexPrism slice 1](../images/primitives3d/sdHexPrism_slice.png) | ![sdHexPrism slice 2](../images/primitives3d/sdHexPrism_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded/fixed by Dillon Cower (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdHexPrism.glsl`](../../glsl/primitives3d/sdHexPrism.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

By the Fundamental Sector Interior Dominance Theorem, inside the $30^\circ$ folded sector of a regular hexagon, $p_y - h_x \le 0 \implies |p_x| - k_z h_x \le 0$ automatically. Thus the 2D edge projection and $z$-cap combine into a single 3D corner `length(max(d, 0.0))`.

---

## GLSL Implementation

```glsl
float sdHexPrism( vec3 p, vec2 h )
{
    const vec3 k = vec3(-0.8660254, 0.5, 0.57735027);
    p = abs(p);
    p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
    vec3 d = vec3(abs(p.x) - k.z*h.x, p.y - h.x, p.z - h.y);
    return min(max(d.y, d.z), 0.0) + length(max(d, 0.0));
}
```

---

## Improvement & Bug-Fix Analysis (Dillon Cower (sdf-gems, 2026))

- **Before (IQ `sdHexPrism`)**: 2 `sqrt`s, 1 `clamp`, 1 `sign`.
- **After (`sdf-gems` `sdHexPrism`)**: **1 `sqrt`, 0 `clamp`s, 0 `sign`s**, `0.0` numerical difference across all of $\mathbb{R}^3$.

### Code Diff vs. Baseline (Inigo Quilez)
```diff
-    vec2 d = vec2(
-       length(p.xy - vec2(clamp(p.x, -k.z*h.x, k.z*h.x), h.x))*sign(p.y - h.x),
-       p.z - h.y );
-    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
+    vec3 d = vec3(abs(p.x) - k.z*h.x, p.y - h.x, p.z - h.y);
+    return min(max(d.y, d.z), 0.0) + length(max(d, 0.0));
```
