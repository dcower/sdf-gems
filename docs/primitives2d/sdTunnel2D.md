# `sdTunnel2D` & `sdgTunnel2D` — Exact 0/1-`sqrt` 2D Arch Tunnel & Analytical Gradient

| 2D SDF Field (Isolines & Zero-Set) |
| :---: |
| ![sdTunnel2D](../images/primitives2d/sdTunnel2D.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions2d/](https://iquilezles.org/articles/distfunctions2d/) (`flSSDy`) |
| **Category** | `primitives2d` |
| **GLSL File** | [`glsl/primitives2d/sdTunnel2D.glsl`](../../glsl/primitives2d/sdTunnel2D.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1` everywhere; `0` `sqrt`s in rectangular strips, `1` `sqrt` in arch & corner; down from `2` `sqrt`s, `4.44e-16` max error) |

---

## Mathematical Overview

### 1. `sdTunnel2D` & `sdgTunnel2D` (`flSSDy`, `2` `sqrt`s $\to$ `0/1` `sqrt`)
Inigo Quilez's `sdTunnel` calls **2 square roots** in the semicircular arch half-plane ($p_y < 0$ after $p_y \gets -p_y$): `length(p)` to compute `d2 = dot2(length(p) - wh.x)` followed by `sqrt(min(d1, d2))`.
- When $p_y > 0$, the semicircular arch ($y \le 0$) is never closer than the rectangular walls, so the signed distance is identical to `sdBox2D`: `(q.x > 0 && q.y > 0) ? length(q) : max(q.x, q.y)` (**0 `sqrt`s** in the rectangular interior and wall strips, **1 `sqrt`** in the bottom corner quadrant).
- When $p_y \le 0$, $q_y = p_y - \text{wh}_y < 0$ always. Outside the arch ($\|p\| \ge \text{wh}_x$), $(\|p\| - \text{wh}_x)^2 \le d_1$ by the triangle inequality; inside the arch ($\|p\| < \text{wh}_x$), $\max(q_x, 0) = 0$ so $-\sqrt{\min(d_1, d_2)} = \max(\|p\| - \text{wh}_x, q_y)$. Thus **both** interior and exterior cases of $p_y \le 0$ collapse to `max(length(p) - wh.x, q.y)`, requiring **only 1 `sqrt`** (`4.44e-16` max error).

![2D Cross & Tunnel Comparison](../images/comparisons/cross2d_tunnel_comparison.png)

---

## GLSL Implementation

```glsl
float sdTunnel2D( vec2 p, vec2 wh )
{
    p.x = abs(p.x); p.y = -p.y;
    vec2 q = p - wh;
    if( p.y > 0.0 )
        return (q.x > 0.0 && q.y > 0.0) ? length(q) : max(q.x, q.y);
    return max(length(p) - wh.x, q.y);
}
```
