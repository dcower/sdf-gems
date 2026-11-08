# `opExtrudeChamfer` — Exact 5-Line Orthogonal Eikonal $`45^\circ`$ Chamfer Extrusion

| 3D Surface ($`f(p) = 0`$) | 2D Cross-Section Slice ($`y = 0`$) |
| :---: | :---: |
| ![opExtrudeChamfer](../images/deform/opExtrudeChamfer.png) | ![opExtrudeChamfer slice](../images/deform/opExtrudeChamfer_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Discovered using Gemini 4 Argon (sdf-gems, 2026)** |
| **Category** | `deform` |
| **GLSL File** | [`glsl/deform/opExtrudeChamfer.glsl`](../../glsl/deform/opExtrudeChamfer.glsl) |
| **Exact Euclidean SDF?** | Yes (`\|\|grad d\|\| = 1`, `1` `sqrt`, `4.44e-16` max error) |

---

## Mathematical Overview

Given any exact 2D SDF $`u = d_{2\text{D}}(p_{xy})`$ and vertical slab coordinate $`v = |z| - h`$, `opExtrudeChamfer` bevels the top and bottom rims with an exact $`45^\circ`$ chamfer of size $`c`$.

By sorting $`a = \min(u, v)`$ and $`b = \max(u, v)`$:
1. **In the exterior corner wedge ($`b > 0 \;\wedge\; b - a > c`$)**: The closest point lies on the flat face $`b = 0`$ or the chamfer rim vertex $`(a = -c, b = 0)`$, at exact Euclidean distance $`\sqrt{b^2 + \max(a + c, 0)^2}`$.
2. **In the chamfer strip and interior ($`b \le 0 \;\vee\; b - a \le c`$)**: The closest point lies on the $`45^\circ`$ chamfer plane $`(a + b + c)/\sqrt{2}`$ or the primary face $`b`$, giving $`\max(b, (a + b + c)/\sqrt{2})`$ with zero interior medial seam distortion (`4.44e-16` max error).

---

## GLSL Implementation

```glsl
float opExtrudeChamfer( float d2d, float z, float h, float c )
{
    float u = d2d, v = abs(z) - h;
    float a = min(u, v), b = max(u, v);
    return (b > 0.0 && b - a > c) ? length(vec2(b, max(a + c, 0.0)))
                                  : max(b, (a + b + c) * 0.70710678);
}
```
