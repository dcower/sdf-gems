# `opUnionStairs` — Stepped Architectural Stairs Union

![opUnionStairs](../images/booleans/opUnionStairs.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Mercury Demogroup** |
| **Source** | [https://mercury.sexy/hg_sdf/](https://mercury.sexy/hg_sdf/) |
| **Category** | `booleans` |
| **GLSL File** | [`glsl/booleans/opUnionStairs.glsl`](../../glsl/booleans/opUnionStairs.glsl) |
| **Exact Euclidean SDF?** | Bound |

---

## Mathematical Overview

Creates $n$ stepped tiers across a blend region of radius $r$ between two SDFs.

---

## GLSL Implementation

```glsl
float opUnionStairs( float a, float b, float r, float n )
{
    float s = r/n;
    float u = b - r;
    return min(min(a, b), 0.5*(u + a + abs((mod(u - a + s, 2.0*s)) - s)));
}
```
