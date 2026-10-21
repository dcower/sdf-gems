# `sdDodecahedron` — 3D Regular Dodecahedron (GDF Plane Bound)

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdDodecahedron](../images/polyhedra3d/sdDodecahedron.png) | ![sdDodecahedron slice 1](../images/polyhedra3d/sdDodecahedron_slice.png) | ![sdDodecahedron slice 2](../images/polyhedra3d/sdDodecahedron_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Mercury Demogroup** |
| **Source** | [https://mercury.sexy/hg_sdf/](https://mercury.sexy/hg_sdf/) |
| **Category** | `polyhedra3d` |
| **GLSL File** | [`glsl/polyhedra3d/sdDodecahedron.glsl`](../../glsl/polyhedra3d/sdDodecahedron.glsl) |
| **Exact Euclidean SDF?** | No (Plane Bound) |

---

## Mathematical Overview

Evaluates the maximum of the 12 pentagonal face planes of a regular dodecahedron.

---

## GLSL Implementation

```glsl
float sdDodecahedron( vec3 p, float r )
{
    p = abs(p);
    const vec2 n = vec2(0.52573111, 0.85065081);
    return max(dot(p.xy, n), max(dot(p.yz, n), dot(p.zx, n))) - r;
}
```
