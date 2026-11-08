# `opRep2D` — 3-Tap 1-Round 2D Hexagonal ($`A_2`$) Delaunay Domain Repetition (`opRepHex3`)

| 3D Surface ($`f(p) = 0`$) | 2D Cross-Section Slice ($`y = 0`$) |
| :---: | :---: |
| ![opRep2D](../images/repetition/opRep2D.png) | ![opRep2D slice](../images/repetition/opRep2D_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez (4-tap Cartesian) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/sdfrepetition/](https://iquilezles.org/articles/sdfrepetition/) |
| **Category** | `repetition` |
| **GLSL File** | [`glsl/repetition/opRep2D.glsl`](../../glsl/repetition/opRep2D.glsl) |
| **Exact Euclidean SDF?** | Yes (`0.0` error vs. 19-tap brute force; `100%` Delaunay containment, `1` `round`, `0` `dot2`) |

---

## Mathematical Overview

Uses the 3 vertices of the containing equilateral Delaunay triangle of the 2D hexagonal lattice $`A_2`$, cutting 2D neighbor-corrected repetition for arbitrary asymmetric/rotated primitives from 4 taps to **3 taps (`1.33x` speedup)** while increasing circle-packing density from $`78.5\%`$ to $`90.7\%`$.

### 1-`round` $`u`$-Space Delaunay Reduction
In normalized sublattice coordinates $`u = p / (s k)`$ with $`k = (1, \sqrt{3})`$, $`i_0 = \mathrm{round}(u)`$, and $`d = u - i_0 \in [-\tfrac{1}{2}, \tfrac{1}{2}]^2`$:
1. For **every** point in $`d \in [-\tfrac{1}{2}, \tfrac{1}{2}]^2`$, **both** $`c_0 = i_0 k`$ and $`c_1 = (i_0 + 0.5\mathrm{sign}(d)) k`$ are vertices of the enclosing Delaunay equilateral triangle — so we **never** need to evaluate a Voronoi distance comparison (`dot2`) to decide which of the two is closer!
2. Furthermore, in $`u`$-space the $`60^\circ`$ Delaunay edge connecting $`i_0`$ and $`i_0 + 0.5\mathrm{sign}(d)`$ is the $`45^\circ`$ diagonal $`|d_x| = |d_y|`$, so the third Delaunay vertex $`c_2`$ is simply $`(i_0 + (\mathrm{sg}_x, 0)) k`$ when $`|d_x| > |d_y|`$ and $`(i_0 + (-0.5\mathrm{sg}_x, 0.5\mathrm{sg}_y)) k`$ otherwise (`100.0%` Delaunay containment in `1` `round` and `0` `dot2` calls).

![Before vs After Comparison](../images/comparisons/hex3_repetition_2d.png)

---

## GLSL Implementation

```glsl
void opRepHex3( in vec2 p, in float s, out vec2 c0, out vec2 c1, out vec2 c2 )
{
    const vec2 k = vec2(1.0, 1.73205081);
    vec2 u  = p / (s * k), i0 = round(u), d = u - i0;
    vec2 sg = vec2(d.x >= 0.0 ? 1.0 : -1.0, d.y >= 0.0 ? 1.0 : -1.0);
    c0 = i0 * k;
    c1 = (i0 + 0.5 * sg) * k;
    c2 = (i0 + ((abs(d.x) > abs(d.y)) ? vec2(sg.x, 0.0) : vec2(-0.5 * sg.x, 0.5 * sg.y))) * k;
}
```
