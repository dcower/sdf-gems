// ============================================================================
// sdHexPrism & sdgHexPrism — Exact 0/1-Sqrt 3D Hexagonal Prism & Analytical Gradient
// @author Inigo Quilez — upgraded to 0/1-sqrt & analytical gradient using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in interior, caps & 6 side prisms, `1` in edge/corner wedges)
// ============================================================================

#ifndef SDF_GEMS_SD_HEX_PRISM_GLSL
#define SDF_GEMS_SD_HEX_PRISM_GLSL

float sdHexPrism( vec3 p, vec2 h )
{
    const vec3 k = vec3(-0.8660254, 0.5, 0.57735027);
    p = abs(p);
    p.xy -= 2.0 * min(dot(k.xy, p.xy), 0.0) * k.xy;
    vec3 d = vec3(abs(p.x) - k.z * h.x, p.y - h.x, p.z - h.y);
    return (max(d.x, min(d.y, d.z)) <= 0.0) ? max(d.y, d.z) : length(max(d, 0.0));
}

#ifndef SDF_GEMS_HAS_SDG_HEX_PRISM
#define SDF_GEMS_HAS_SDG_HEX_PRISM
vec4 sdgHexPrism( vec3 p, vec2 h )
{
    const vec3 k = vec3(-0.86602540, 0.5, 0.57735027);
    vec3 s = vec3((p.x < 0.0) ? -1.0 : 1.0,
                  (p.y < 0.0) ? -1.0 : 1.0,
                  (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    float dk = dot(k.xy, p.xy);
    bool  rf = (dk < 0.0);
    if( rf ) p.xy -= 2.0 * dk * k.xy;

    vec3  d = vec3(abs(p.x) - k.z * h.x, p.y - h.x, p.z - h.y);
    float dist;
    vec3  g;
    if( max(d.x, min(d.y, d.z)) <= 0.0 )
    {
        dist = max(d.y, d.z);
        g    = (d.z > d.y) ? vec3(0.0, 0.0, 1.0) : vec3(0.0, 1.0, 0.0);
    }
    else
    {
        vec3  v    = max(d, 0.0);
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        dist = v2 * invV;
        g    = vec3((p.x < 0.0) ? -v.x : v.x, v.y, v.z) * invV;
    }
    if( rf ) g.xy -= 2.0 * dot(k.xy, g.xy) * k.xy;
    return vec4(dist, s * g);
}
#endif // SDF_GEMS_HAS_SDG_HEX_PRISM

#endif // SDF_GEMS_SD_HEX_PRISM_GLSL
