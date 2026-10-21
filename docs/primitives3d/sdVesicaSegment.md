# `sdVesicaSegment` — Exact 3D Revolved Vesica / Spindle

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdVesicaSegment](../images/primitives3d/sdVesicaSegment.png) | ![sdVesicaSegment slice](../images/primitives3d/sdVesicaSegment_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdVesicaSegment.glsl`](../../glsl/primitives3d/sdVesicaSegment.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Revolves a 2D vesica around its chord $[a, b]$ to form a 3D spindle.

---

## GLSL Implementation

```glsl
float sdVesicaSegment( vec3 p, vec3 a, vec3 b, float w )
{
    vec3 c = (a + b)*0.5;
    float l = length(b - a);
    vec3 v = (b - a)/l;
    float y = dot(p - c, v);
    vec2 q = vec2(length(p - c - y*v), abs(y));
    float r = 0.5*l;
    float d = 0.5*(r*r - w*w)/w;
    vec3 h = (r*q.x < d*(q.y - r)) ? vec3(0.0, r, 0.0) : vec3(-d, 0.0, d + w);
    return length(q - h.xy) - h.z;
}
```
