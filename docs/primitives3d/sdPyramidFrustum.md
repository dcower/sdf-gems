# `sdPyramidFrustum` — Exact 1-Sqrt 3D Square Frustum

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdPyramidFrustum](../images/primitives3d/sdPyramidFrustum.png) | ![sdPyramidFrustum slice 1](../images/primitives3d/sdPyramidFrustum_slice.png) | ![sdPyramidFrustum slice 2](../images/primitives3d/sdPyramidFrustum_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Dillon Cower (sdf-gems, 2026)** |
| **Source** | [https://github.com/dcower/sdf-gems](https://github.com/dcower/sdf-gems) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdPyramidFrustum.glsl`](../../glsl/primitives3d/sdPyramidFrustum.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Computes the exact 3D Euclidean SDF of a square frustum of bottom half-width $b_1$, top half-width $b_2$, and half-height $h$ in a **single `sqrt`**.

---

## GLSL Implementation

```glsl
float sdPyramidFrustum( vec3 p, float b1, float b2, float h )
{
    p.xz = abs(p.xz);
    p.xz = (p.z > p.x) ? p.zx : p.xz;
    float bm = 0.5*(b1 + b2), k = 0.5*(b1 - b2);
    vec3  q  = p - vec3(bm, 0.0, bm);
    float u1 = clamp((h*q.y - k*q.x)/(h*h + k*k), -1.0, 1.0);
    float u2 = clamp((h*q.y - k*(q.x + q.z))/(h*h + 2.0*k*k), -1.0, 1.0);
    float d1 = dot2(vec3(q.x + k*u1, q.y - h*u1, max(q.z + k*u1, 0.0)));
    float d2 = dot2(vec3(q.x + k*u2, q.y - h*u2, q.z + k*u2));
    float d  = min(d1, d2);
    if( p.x < ((p.y < 0.0) ? b1 : b2) ) d = min(d, (abs(p.y) - h)*(abs(p.y) - h));
    return sqrt(d) * ((max(h*q.x + k*q.y, abs(q.y) - h) < 0.0) ? -1.0 : 1.0);
}
```
