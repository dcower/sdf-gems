# `sdDeathStar` — Exact 2-Length 3D Sphere Subtracted by Sphere (Confocal Ellipse Bisector Reduction)

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdDeathStar](../images/primitives3d/sdDeathStar.png) | ![sdDeathStar slice](../images/primitives3d/sdDeathStar_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdDeathStar.glsl`](../../glsl/primitives3d/sdDeathStar.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, `2` dynamic `length()` calls down from `3`, `0.0` error) |

---

## Mathematical Overview & Confocal Conic Bisector Theorem

Inigo Quilez's baseline `sdDeathStar` (and 2D `sdMoon`) evaluates **two separate 2D `length()` calls** in its `else` branch:
$$\max\big(\|p\| - r_a, \; -( \|p - (d, 0)\| - r_b )\big)$$
because for unequal radii $r_a \ne r_b$ (or shifted centers), the medial axis where the signed distances to the two spheres are equal is a curved 4th-degree algebraic surface rather than a flat coordinate plane.

We prove the **Confocal Conic Bisector Theorem**, which decides which sphere is closer using a **sqrt-free quadratic test** and collapses all three Voronoi regions (the circular rim, Sphere $A$, and Crater Sphere $B$) into a **single `length(p - c.xy)` call**:

1. **Confocal Ellipse Medial Axis**:
   Let $r_1 = \|p\|$ and $r_2 = \|p - (d, 0)\|$. The outer sphere $S_A(0, r_a)$ dominates the crater sphere $S_B((d, 0), r_b)$ iff:
   $$r_1 - r_a \ge r_b - r_2 \iff r_1 + r_2 \ge R, \qquad R = r_a + r_b > 0.$$
   Geometrically, the locus $r_1 + r_2 = R$ is a **confocal ellipse** with foci at the two sphere centers $(0, 0)$ and $(d, 0)$ that passes directly through the rim intersection points $(a, \pm b)$.
2. **Exact Radial Cancellation ($r_2^2 - r_1^2 = d^2 - 2 d p_x$)**:
   Because both spheres share $p_y^2$, the difference of squared radii $r_2^2 - r_1^2 = (p_x - d)^2 - p_x^2 = d^2 - 2 d p_x$ is **linear in $p_x$**. Using the rim abscissa $a = \frac{r_a^2 - r_b^2 + d^2}{2d}$ (so $d^2 = 2 d a - r_a^2 + r_b^2$), squaring $r_2 \ge R - r_1$ cancels $r_1^2$ identically:
   $$r_2^2 - (R - r_1)^2 = (d^2 - 2 d p_x) - R^2 + 2 R r_1 = 2\big(R r_1 - \underbrace{(r_a R + d(p_x - a))}_{m_0}\big).$$
3. **Universal Clamped Equivalence**:
   When $r_1 \ge R$ ($R r_1 \ge R^2$), $r_1 + r_2 \ge R$ holds trivially since $r_2 \ge 0$; when $r_1 < R$, both sides of $r_2 \ge R - r_1$ are non-negative so squaring is an exact equivalence. Defining $m = \operatorname{clamp}(r_a R + d(p_x - a), \, 0, \, R^2)$ yields the universal sqrt-free condition across all $(r_a, r_b, d)$:
   $$\boxed{\|p\| - r_a \ge r_b - \|p - (d, 0)\| \iff R^2 \|p\|^2 \ge m^2.}$$

---

## GLSL Implementation

```glsl
float sdDeathStar( vec3 p2, float ra, float rb, float d )
{
    float a = (ra*ra - rb*rb + d*d)/(2.0*d);
    float b = sqrt(max(ra*ra - a*a, 0.0));
    vec2  p = vec2( p2.x, length(p2.yz) );
    float r = ra + rb, m = clamp(ra*r + d*(p.x - a), 0.0, r*r);
    vec3  c = (p.x*b - p.y*a > d*max(b - p.y, 0.0)) ? vec3(a, b, 0.0) :
              (r*r*dot(p, p) >= m*m) ? vec3(0.0, 0.0, ra) : vec3(d, 0.0, -rb);
    return length(p - c.xy) * ((c.z < 0.0) ? -1.0 : 1.0) - c.z;
}
```

---

## Improvement & Optimization Analysis (Discovered using Gemini 4 Argon (sdf-gems, 2026))

![Confocal Conic Bisector Reduction](../images/comparisons/death_star_comparison.png)

| Metric | Baseline `sdDeathStar` | `sdf-gems` `sdDeathStar` |
| :--- | :---: | :---: |
| **Dynamic `length()` (`sqrt`) Calls (3D)** | `3` | **`2` (`33%` fewer dynamic `sqrt`s)** |
| **Dynamic `length()` (`sqrt`) Calls (2D `sdMoon2D`)** | `2` | **`1` (`50%` fewer dynamic `sqrt`s)** |
| **Max Error vs. Baseline** | — | **`0.0` (`0.000e+00` in FP64 across $5 \times 10^6$ points)** |

### Code Diff vs. Baseline (Inigo Quilez)
```diff
 float sdDeathStar( vec3 p2, float ra, float rb, float d )
 {
     float a = (ra*ra - rb*rb + d*d)/(2.0*d);
     float b = sqrt(max(ra*ra - a*a, 0.0));
     vec2 p = vec2( p2.x, length(p2.yz) );
-    if( p.x*b - p.y*a > d*max(b - p.y, 0.0) )
-        return length(p - vec2(a,b));
-    else
-        return max( (length(p) - ra), -(length(p - vec2(d,0.0)) - rb));
+    float r = ra + rb, m = clamp(ra*r + d*(p.x - a), 0.0, r*r);
+    vec3  c = (p.x*b - p.y*a > d*max(b - p.y, 0.0)) ? vec3(a, b, 0.0) :
+              (r*r*dot(p, p) >= m*m) ? vec3(0.0, 0.0, ra) : vec3(d, 0.0, -rb);
+    return length(p - c.xy) * ((c.z < 0.0) ? -1.0 : 1.0) - c.z;
 }
```
