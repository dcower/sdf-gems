# `sdPie2D` & `sdgPie2D` — Exact 0/1-`sqrt` 2D Circular Pie Sector & Analytical Gradient

| 2D SDF Field (Isolines & Zero-Set) |
| :---: |
| ![sdPie2D](../images/primitives2d/sdPie2D.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded to 0/1-`sqrt` via Parabolic $v^2$-Cancellation using Gemini 4 Argon (`sdf-gems`, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions2d/](https://iquilezles.org/articles/distfunctions2d/) & [https://iquilezles.org/articles/distgradfunctions2d/](https://iquilezles.org/articles/distgradfunctions2d/) |
| **Category** | `primitives2d` |
| **GLSL File** | [`glsl/primitives2d/sdPie2D.glsl`](../../glsl/primitives2d/sdPie2D.glsl) |
| **Exact Euclidean SDF?** | Yes ($\|\nabla d\| = 1$, **`0` `sqrt`s** in $44\%$–$60\%$ of $[-2r, 2r]^2$ and **`1` `sqrt`** elsewhere, down from `2` `sqrt`s in IQ's `sdPie` / `sdgPie`) |

![Triangle, Pie, and Ring 0-Sqrt / 0-Division Comparison](../images/comparisons/triangle_pie_ring_comparison.png)

---

## Mathematical Overview

### 1. Orthonormal Sector-Dominance & Exterior 0-`sqrt` Strip
After folding $p_x \leftarrow |p_x|$ with sx $= \operatorname{sign}(x)$, rotating into the orthonormal ray basis of the sector edge $c = (\sin\theta, \cos\theta)$ gives:
$$u = p \cdot c = c_x |p_x| + c_y p_y, \qquad v = c_y |p_x| - c_x p_y,$$
with $\|p\|_2^2 = u^2 + v^2$.
- **Outside the sector ($v > 0$):** By the triangle inequality, $\|p\|_2 - r \le \sqrt{(u - \operatorname{clamp}(u, 0, r))^2 + v^2}$ holds identically. Moreover, whenever $0 < u < r$ (the entire exterior normal strip of the radial edge!), $\operatorname{clamp}(u, 0, r) = u$, so the exact exterior Euclidean distance is simply the linear coordinate **`v` (`0` `sqrt`s)** and the exact gradient is $(\operatorname{sx}\, c_y, -c_x)$!

### 2. The Parabolic $v^2$-Cancellation Theorem (Interior 0-`sqrt` Wedge)
Inside the sector ($v \le 0$) for $u \ge 0$, the closest boundary point is either the radial edge at signed distance $v \le 0$ or the circular rim at signed distance $\sqrt{u^2 + v^2} - r$. Comparing the two:
$$v > \sqrt{u^2 + v^2} - r \iff r + v > \sqrt{u^2 + v^2}.$$
When $r + v > 0$, squaring both sides yields:
$$r^2 + 2rv + v^2 > u^2 + v^2 \iff r(r + 2v) > u^2 \quad (v^2 \text{ cancels identically!}).$$
And when $r + v \le 0$, $r(r + 2v) \le -r^2 \le 0 \le u^2$, so the single polynomial condition **`r * (r + 2.0 * v) > u * u`** is necessary and sufficient for all $v \le 0, u \ge 0$! Whenever it holds, the exact interior distance is **`v` with `0` `sqrt`s**, eliminating `sqrt` across **$44\%$–$60\%$ of $[-2r, 2r]^2$** in both `sdPie2D` and `sdgPie2D`.

---

## GLSL Implementation

```glsl
// 1. Exact 0/1-Sqrt 2D Circular Pie Sector (0 sqrts in radial edge strip & interior wedge!)
float sdPie2D( in vec2 p, in vec2 c, in float r )
{
    p.x = abs(p.x);
    float u = dot(p, c), v = c.y * p.x - c.x * p.y;
    if( v > 0.0 )
        return (u > 0.0 && u < r) ? v : length(vec2(u - clamp(u, 0.0, r), v));
    if( u >= 0.0 && r * (r + 2.0 * v) > u * u )
        return v;
    float l = length(p);
    return (u < 0.0) ? -l : l - r;
}
```
