# `opExtrudeMolding` — Exact 1-Sqrt Cove, Bead & Roman Ogee Molding Extrusions

![opExtrudeMolding](../images/deform/opExtrudeMolding.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://github.com/dcower/sdf-gems](https://github.com/dcower/sdf-gems) |
| **Category** | `deform` |
| **GLSL File** | [`glsl/deform/opExtrudeMolding.glsl`](../../glsl/deform/opExtrudeMolding.glsl) |
| **Exact Euclidean SDF?** | Yes (`\|\|grad d\|\| = 1`, `1` `sqrt` total, `0.0` to `4.44e-16` max error) |

---

## Mathematical Overview

![Exact 1-Sqrt Architectural Molding Extrusions](../images/comparisons/extrude_molding_comparison.png)

By the **Orthogonal Eikonal Composition Theorem**, given $`(u, v) = (d_{\text{2D}}(p_{xy}),\; |z| - h)`$ and sorted coordinates $`(a, b) = (\min(u, v), \max(u, v))`$, any exact 2D corner profile lifts to an exact 3D molding extrusion with $`\|\nabla_{\text{3D}} F\| \equiv 1`$:
1. **`opExtrudeCove` (Concave Quarter-Circle Rabbet / Cavetto)**:
   Carves a quarter-circle concave cove of radius $`r`$ centered at $`(0, 0)`$ out of the prism corner $`u \le 0, v \le 0`$. For $`b > 0`$, the closest point on the third-quadrant arc is always its endpoint $`(-r, 0)`$, giving $`d = \sqrt{\max(a + r, 0)^2 + b^2}`$; for $`b \le 0`$, $`d = \max(b, r - \sqrt{a^2 + b^2})`$ (**4 lines, 1 `sqrt`**).
2. **`opExtrudeBead` (Convex Astragal Torus Bead — Reflection Dual of `opExtrudeCove`)**:
   Unions a radius-$`r`$ circular bead centered at $`(0, 0)`$ with $`\{u \le 0, v \le 0\}`$, eliminating the $`29.3\%`$ interior quirk-corner error of naive `min(d_extrude, length(vec2(u, v)) - r)` in **4 lines and 1 `sqrt`**.
3. **`opExtrudeOgee` (Tangent Roman Ogee / Cyma Recta)**:
   Spans $`[-r, 0]^2`$ with an upper concave quarter-circle $`A_1`$ centered at $`(-R, 0)`$ and a lower convex quarter-circle $`A_2`$ centered at $`(-R, -r)`$ ($`R = r/2`$), meeting with $`C^1`$ horizontal tangency at $`M = (-R, -R)`$. Squaring the Apollonius equidistant equation $`\sqrt{x^2 + (y + R)^2} - \sqrt{(x + R)^2 + (y - R)^2} = R`$ cancels $`(x + R)^2`$ identically into the pure quadratic condition $`3 y^2 - 2 y (2 x + R) > R^2`$, selecting the active center **before** taking the single `sqrt`.

| Molding Operator | `sqrt` Count | Max SDF Error vs Ground Truth | Max Eikonal Error $`\big\vert \Vert \nabla d\Vert - 1\big\vert`$ | Naive Boolean Error |
| :--- | :---: | :---: | :---: | :---: |
| `opExtrudeCove` / `sdgExtrudeCove` | $`1`$ | `0.000e+00` | `2.220e-16` | `0.3142` |
| `opExtrudeBead` / `sdgExtrudeBead` | $`1`$ | `0.000e+00` | `2.220e-16` | `0.1050` ($`29.3\%`$ interior) |
| `opExtrudeOgee` / `sdgExtrudeOgee` | $`0\text{-}1`$ | `4.441e-16` | `2.220e-16` | `0.2828` |

---

## GLSL Implementation

```glsl
float opExtrudeCove( in float d2d, in float z, in float h, in float r )
{
    float u = d2d, v = abs(z) - h;
    float a = min(u, v), b = max(u, v);
    float l = length((b > 0.0) ? vec2(max(a + r, 0.0), b) : vec2(a, b));
    return (b > 0.0) ? l : max(b, r - l);
}

float opExtrudeBead( in float d2d, in float z, in float h, in float r )
{
    float u = d2d, v = abs(z) - h;
    float a = min(u, v), b = max(u, v);
    float l = length((b < 0.0) ? vec2(max(a + r, 0.0), b) : vec2(a, b));
    return (b < 0.0) ? -l : ((a > -r) ? l - r : min(b, l - r));
}

float opExtrudeOgee( in float d2d, in float z, in float h, in float r )
{
    float u = d2d, v = abs(z) - h, R = 0.5 * r;
    float x = u + R, y = v + R;
    if( x > 0.0 && y <= -R ) return u;
    bool  c = (x <= 0.0) ? (y > R) : (y > 0.0 && 3.0 * y * y - 2.0 * y * (2.0 * x + R) > R * R);
    vec2  q = vec2(c ? max(x + R, 0.0) : x, (c || x <= 0.0) ? y - R : y + R);
    float l = length(q);
    return c ? l : ((x <= 0.0) ? max(max(u, v), R - l) : l - R);
}
```
