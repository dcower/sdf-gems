# `sdSegment2D` — Exact 2D Line Segment / Capsule

![sdSegment2D](../images/primitives2d/sdSegment2D.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions2d/](https://iquilezles.org/articles/distfunctions2d/) |
| **Category** | `primitives2d` |
| **GLSL File** | [`glsl/primitives2d/sdSegment2D.glsl`](../../glsl/primitives2d/sdSegment2D.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Exact Euclidean distance to a 2D line segment $[a, b]$.

---

## GLSL Implementation

```glsl
float sdSegment2D( vec2 p, vec2 a, vec2 b )
{
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}
```
