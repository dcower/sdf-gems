# `sdOctahedron` — Exact 3D Regular Octahedron

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdOctahedron](../images/polyhedra3d/sdOctahedron.png) | ![sdOctahedron slice 1](../images/polyhedra3d/sdOctahedron_slice.png) | ![sdOctahedron slice 2](../images/polyhedra3d/sdOctahedron_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `polyhedra3d` |
| **GLSL File** | [`glsl/polyhedra3d/sdOctahedron.glsl`](../../glsl/polyhedra3d/sdOctahedron.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Computes the exact Euclidean distance to a regular octahedron of half-span $s$ via octant folding and 3-way edge selection.

---

## GLSL Implementation

```glsl
float sdOctahedron( vec3 p, float s )
{
    p = abs(p);
    float m = p.x + p.y + p.z - s;
    vec3 q;
         if( 3.0*p.x < m ) q = p.xyz;
    else if( 3.0*p.y < m ) q = p.yzx;
    else if( 3.0*p.z < m ) q = p.zxy;
    else return m*0.57735027;

    float k = clamp(0.5*(q.z - q.y + s), 0.0, s);
    return length(vec3(q.x, q.y - s + k, q.z - k));
}
```
