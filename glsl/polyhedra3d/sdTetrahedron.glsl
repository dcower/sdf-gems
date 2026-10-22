// ============================================================================
// sdTetrahedron & sdgTetrahedron — Exact 0/1-Sqrt 3D Regular Tetrahedron (Td Parity-Sorted Weyl Chamber)
// @author Mercury Demogroup (plane bound) — upgraded to 0/1-sqrt & analytical gradient using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in 100% of interior & 4 face prisms, `1` sqrt in ridge/vertex wedges)
// ============================================================================

#ifndef SDF_GEMS_SD_TETRAHEDRON_GLSL
#define SDF_GEMS_SD_TETRAHEDRON_GLSL

float sdTetrahedron( vec3 p, float s )
{
    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;
    p.z *= sig;

    float m = p.x + p.y - p.z - s;
    if( m <= 1.5 * (p.y - p.z) )
        return m * 0.57735027; // 0 sqrts in 100% of interior & all 4 exterior face prisms!
    float t = min(0.5 * (p.y + p.z), s);
    return length(vec3(p.x - s, p.y - t, p.z - t));
}

vec4 sdgTetrahedron( vec3 p, float s )
{
    vec3 sg = vec3((p.x < 0.0) ? -1.0 : 1.0,
                   (p.y < 0.0) ? -1.0 : 1.0,
                   (p.z < 0.0) ? -1.0 : 1.0);
    float sig = sg.x * sg.y * sg.z;
    p = abs(p);
    bool s1 = (p.x < p.y); if( s1 ) p.xy = p.yx;
    bool s2 = (p.y < p.z); if( s2 ) p.yz = p.zy;
    bool s3 = (p.x < p.y); if( s3 ) p.xy = p.yx;
    p.z *= sig;

    float m = p.x + p.y - p.z - s;
    float d;
    vec3  g;
    if( m <= 1.5 * (p.y - p.z) )
    {
        d = m * 0.57735027;
        g = vec3(0.57735027, 0.57735027, -0.57735027);
    }
    else
    {
        float t    = min(0.5 * (p.y + p.z), s);
        vec3  v    = vec3(p.x - s, p.y - t, p.z - t);
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = v * invV;
    }

    g.z *= sig;
    if( s3 ) g.xy = g.yx;
    if( s2 ) g.yz = g.zy;
    if( s1 ) g.xy = g.yx;
    return vec4(d, sg * g);
}

#endif // SDF_GEMS_SD_TETRAHEDRON_GLSL
