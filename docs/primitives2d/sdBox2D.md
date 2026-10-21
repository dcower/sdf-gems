# `sdBox2D` — Exact 2D Box

![sdBox2D](../images/primitives2d/sdBox2D.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions2d/](https://iquilezles.org/articles/distfunctions2d/) |
| **Category** | `primitives2d` |
| **GLSL File** | [`glsl/primitives2d/sdBox2D.glsl`](../../glsl/primitives2d/sdBox2D.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Exact Euclidean signed distance to a 2D rectangle of half-extents $b$.

---

## GLSL Implementation

```glsl
float sdBox2D( vec2 p, vec2 b )
{
    vec2 d = abs(p) - b;
    return length(max(d, 0.0)) + min(max(d.x, d.y), 0.0);
}
```
