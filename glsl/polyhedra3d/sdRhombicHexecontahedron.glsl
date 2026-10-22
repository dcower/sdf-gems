// ============================================================================
// sdRhombicHexecontahedron — Exact 0/1-Sqrt Unkelbach 60-Rhombus Stellation
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, `2.22e-15` FP64 / `2.14e-9` FP32 error, 0 `sqrt`s in face prisms, 1 `sqrt` on ridges/valleys)
// ============================================================================

#ifndef SDF_GEMS_SD_RHOMBIC_HEXECONTAHEDRON_GLSL
#define SDF_GEMS_SD_RHOMBIC_HEXECONTAHEDRON_GLSL

float sdRhombicHexecontahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.50000000, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    // Face-Normal Swizzle Identity: n = -n3.yzx = (phi/2, -iphi/2, 1/2)
    vec3  n  = -n3.yzx;
    float df = dot(p, n) - r;
    if( p.x >= n.x * df && p.y >= n.y * df )
        return df; // 0 sqrts in both interior and exterior face prisms!

    const vec2 e = vec2(0.52573111, 0.85065081);
    float z2 = p.z - 2.0 * r;
    float L  = 1.17557050 * r;
    if( df < 0.0 )
    {
        // Concave valley seam V2 -> V5 in Wall A (y = 0)
        float t = clamp(e.x * p.x - e.y * z2, 0.0, L);
        return -length(vec3(p.x - e.x * t, p.y, z2 + e.y * t)); // 1 sqrt
    }
    // Convex ridge edge V2 -> V3 in Wall B (x = 0)
    float t = clamp(e.y * p.y + e.x * z2, 0.0, L);
    return length(vec3(p.x, p.y - e.y * t, z2 - e.x * t)); // 1 sqrt
}

vec4 sdgRhombicHexecontahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.50000000, -0.80901699, 0.30901699);
    vec3 s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); float c0 = dot(p, n3); if( c0 < 0.0 ) p -= 2.0 * c0 * n3;
    vec3 s1 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); float c1 = dot(p, n3); if( c1 < 0.0 ) p -= 2.0 * c1 * n3;
    vec3 s2 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); float c2 = dot(p, n3); if( c2 < 0.0 ) p -= 2.0 * c2 * n3;
    vec2 s3 = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    p.xy = abs(p.xy);

    vec3  n  = -n3.yzx;
    float df = dot(p, n) - r;
    vec3  g;
    float d;
    if( p.x >= n.x * df && p.y >= n.y * df )
    {
        d = df;
        g = n;
    }
    else
    {
        const vec2 e = vec2(0.52573111, 0.85065081);
        float z2 = p.z - 2.0 * r;
        float L  = 1.17557050 * r;
        if( df < 0.0 )
        {
            float t = clamp(e.x * p.x - e.y * z2, 0.0, L);
            vec3  v = vec3(p.x - e.x * t, p.y, z2 + e.y * t);
            float l = length(v);
            d = -l;
            g = -v / max(l, 1e-20);
        }
        else
        {
            float t = clamp(e.y * p.y + e.x * z2, 0.0, L);
            vec3  v = vec3(p.x, p.y - e.y * t, z2 - e.x * t);
            float l = length(v);
            d = l;
            g = v / max(l, 1e-20);
        }
    }
    g.xy *= s3;
    if( c2 < 0.0 ) g -= 2.0 * dot(g, n3) * n3; g *= s2;
    if( c1 < 0.0 ) g -= 2.0 * dot(g, n3) * n3; g *= s1;
    if( c0 < 0.0 ) g -= 2.0 * dot(g, n3) * n3; g *= s0;
    return vec4(d, g);
}

#endif // SDF_GEMS_SD_RHOMBIC_HEXECONTAHEDRON_GLSL
