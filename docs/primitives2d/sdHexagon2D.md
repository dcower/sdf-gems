# `sdHexagon2D` — Exact 0/1-`sqrt` 2D Regular Hexagon

![sdHexagon2D](../images/primitives2d/sdHexagon2D.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded to 0-`sqrt` interior/face-strip using Gemini 4 Argon (`sdf-gems`, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions2d/](https://iquilezles.org/articles/distfunctions2d/) |
| **Category** | `primitives2d` |
| **GLSL File** | [`glsl/primitives2d/sdHexagon2D.glsl`](../../glsl/primitives2d/sdHexagon2D.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, `0` `sqrt`s in interior & face strips, `1` in corner wedges) |

---

## Mathematical Overview

Folds $\mathbb{R}^2$ into the $30^\circ$ fundamental sector of $D_6$ where $|p_x| \le k_z p_y$ ($k_z = \tan(\pi/6) = 1/\sqrt{3}$). By **Fundamental-Chamber Vertex Dominance**, for $d = (|p_x| - k_z r, p_y - r)$, $d_x > 0 \implies |p_x| > k_z r \implies k_z p_y \ge |p_x| > k_z r \implies d_y > 0$ automatically. Thus `length(d)` is only evaluated in the $6$ exterior corner wedges ($d_x > 0$), while the entire interior and all $6$ flat-edge normal strips ($d_x \le 0$) return $d_y = p_y - r$ directly with **`0` `sqrt`s, `0` `clamp`s, and `0` `sign`s**.

---

## GLSL Implementation

```glsl
float sdHexagon2D( vec2 p, float r )
{
    const vec3 k = vec3(-0.866025404, 0.5, 0.577350269);
    p = abs(p);
    p -= 2.0*min(dot(k.xy, p), 0.0)*k.xy;
    vec2 d = vec2(abs(p.x) - k.z*r, p.y - r);
    return (d.x > 0.0) ? length(d) : d.y;
}
```
