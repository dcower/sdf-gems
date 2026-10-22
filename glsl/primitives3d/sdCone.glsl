// ============================================================================
// sdCone — Exact 0/1/2-Sqrt 3D Finite & 1-Sqrt Infinite Cone (Slanted-Wall Cap-Voronoi Theorem)
// @author Inigo Quilez — upgraded to 0/1/2-sqrt using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in base-cap Voronoi, `1` in apex/mantle, `2` in rim wedge)
// ============================================================================

#ifndef SDF_GEMS_SD_CONE_GLSL
#define SDF_GEMS_SD_CONE_GLSL

// 1. Exact 0/1/2-Sqrt Finite 3D Cone (c = vec2(sin(alpha), cos(alpha)), h = height)
//    By the Slanted-Wall Cap-Voronoi Radius Theorem:
//      rc = (h * c.x + (1.0 + c.x) * min(d, 0.0)) / c.y,   d = -p.y - h,
//    whenever rc >= 0 and dot(p.xz, p.xz) <= rc^2 (covering the entire exterior base-cap
//    cylinder AND the interior base-cap Voronoi frustum!), the exact 3D Euclidean distance
//    is simply `d = -p.y - h` with ZERO square roots!
float sdCone( vec3 p, vec2 c, float h )
{
    float x2 = dot(p.xz, p.xz), d = -p.y - h;
    float rc = (h * c.x + (1.0 + c.x) * min(d, 0.0)) / c.y;
    if( rc >= 0.0 && x2 <= rc * rc )
        return d;
    if( p.y >= 0.0 && c.x * c.x * x2 <= c.y * c.y * p.y * p.y )
        return sqrt(x2 + p.y * p.y);
    float wx = sqrt(x2), rx = wx - h * c.x / c.y;
    return (rx > 0.0 && c.x * rx + c.y * d >= 0.0)
         ? length(vec2(rx, d))
         : c.y * wx + c.x * p.y;
}

// 2. Infinite 3D Cone (1 `sqrt` everywhere in R^3, down from 2 `sqrt`s)
float sdConeInf( vec3 p, vec2 c )
{
    float x2 = dot(p.xz, p.xz);
    return (p.y >= 0.0 && c.x * c.x * x2 <= c.y * c.y * p.y * p.y)
         ? sqrt(x2 + p.y * p.y)
         : sqrt(x2) * c.y + p.y * c.x;
}

#endif // SDF_GEMS_SD_CONE_GLSL
