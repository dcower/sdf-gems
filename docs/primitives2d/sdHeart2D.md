# `sdHeart2D` & `sdgHeart2D` — 0/1-`sqrt` Exact 2D Heart & `NaN`-Free Analytical Gradient

| 2D SDF Field (`sdHeart2D`) |
| :---: |
| ![sdHeart2D](../images/primitives2d/sdHeart2D.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded to 0/1-`sqrt` & `NaN`-free gradient using Gemini 4 Argon (`sdf-gems`, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions2d/](https://iquilezles.org/articles/distfunctions2d/) & [https://iquilezles.org/articles/distgradfunctions2d/](https://iquilezles.org/articles/distgradfunctions2d/) |
| **Category** | `primitives2d` |
| **GLSL File** | [`glsl/primitives2d/sdHeart2D.glsl`](../../glsl/primitives2d/sdHeart2D.glsl) |
| **Exact Euclidean SDF?** | Yes ($`\Vert \nabla d\Vert = 1`$, `0` `sqrt`s in $`24.4\%`$ of $`[-2,2]^2`$, `0` `NaN`s on the $`45^\circ`$ zero-set) |

---

## Mathematical Overview

### 0/1-`sqrt` `sdHeart2D` & `NaN`-Free `sdgHeart2D`
In IQ's `sdHeart`, for $`s = |p_x| + p_y \le 1`$, IQ compares `dot2(p - vec2(0, 1))` against `dot2(p - 0.5*max(s, 0.0))` and takes a `sqrt`. Expanding their difference gives the exact algebraic identity:
```math
\|p - (0,1)\|_2^2 - \frac{1}{2}(p_x - p_y)^2 \equiv \frac{1}{2}(s - 1)^2 + \frac{1}{2} + (p_x - p_y).
```
Whenever $`0 < s \le 1`$ and $`\frac{1}{2}(s - 1)^2 + \frac{1}{2} + (p_x - p_y) \ge 0`$ (which covers **100% of the lower exterior** $`p_x \ge p_y`$ and the interior strip $`p_x - p_y \ge -0.5`$, or **$`24.4\%`$ of $`[-2, 2]^2`$**), the exact signed distance is simply $`(p_x - p_y)\sqrt{0.5}`$ with **0 `sqrt`s** ($`4.44 \times 10^{-16}`$ max error)! This also eliminates the **`0.0 / 0.0 = NaN`** division on the $`45^\circ`$ zero-set $`p_x = p_y \in [0, 0.5]`$ in IQ's `sdgHeart` (`DldXRf`).

*(For `sdEgg2D` / `sdEgg3D` and `sdHorseshoe2D` / `sdHorseshoe3D`, see [`sdEgg2D.md`](sdEgg2D.md) and [`sdHorseshoe2D.md`](sdHorseshoe2D.md).)*

---

## GLSL Implementation

```glsl
// 0/1-Sqrt Exact 2D Heart
float sdHeart2D( in vec2 p )
{
    p.x = abs(p.x);
    float s = p.x + p.y, d = p.x - p.y;
    if( s > 1.0 ) return length(p - vec2(0.25, 0.75)) - 0.35355339;
    if( s > 0.0 && 0.5 * (s - 1.0) * (s - 1.0) + 0.5 + d >= 0.0 ) return d * 0.70710678;
    return (s <= 0.0) ? length(p) : -length(p - vec2(0.0, 1.0));
}
```
