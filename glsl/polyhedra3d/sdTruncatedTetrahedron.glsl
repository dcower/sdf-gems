// ============================================================================
// sdTruncatedTetrahedron & sdgTruncatedTetrahedron — Exact 0/1-Sqrt 3D Archimedean Truncated Tetrahedron
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & face prisms, 1 sqrt in edge/vertex cone)
// ============================================================================

// 1. Exact 0/1-Sqrt Archimedean Truncated Tetrahedron (4 regular hexagons + 4 equilateral triangles):
//    Folds R^3 into the order-24 T_d Weyl chamber x >= y >= |z| using 3 compare-swaps + parity sign.
float sdTruncatedTetrahedron( vec3 p, float s )
{
    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;
    p.z *= sig;

    float c     = (1.0 / 3.0) * s;
    float d_hex = (p.x + p.y - p.z - s) * 0.57735027;
    float d_tri = (p.x + p.y + p.z - 5.0 * c) * 0.57735027;
    if( max(d_hex, d_tri) <= 0.0 || min(p.x + p.y - 2.0 * p.z - 2.0 * c, p.x + p.z + max(p.x - p.y, p.y + p.z) - 2.0 * s) <= 0.0 )
        return max(d_hex, d_tri);

    float t = max(c - 0.5 * min(p.x - p.y, p.y + p.z), 0.0);
    vec3  d = (p.x - p.y <= p.y + p.z)
            ? vec3(p.x - s + t, p.y - c - t, p.z - c)
            : vec3(p.x - s, p.y - c + t, p.z - c + t);
    return length(d);
}

// 2. Analytical 3D Distance and Unit Gradient vec4(d, grad.xyz) (0/1-inversesqrt)
vec4 sdgTruncatedTetrahedron( vec3 p, float s )
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

    float c     = (1.0 / 3.0) * s;
    float d_hex = (p.x + p.y - p.z - s) * 0.57735027;
    float d_tri = (p.x + p.y + p.z - 5.0 * c) * 0.57735027;
    float d;
    vec3  g;
    if( max(d_hex, d_tri) <= 0.0 || min(p.x + p.y - 2.0 * p.z - 2.0 * c, p.x + p.z + max(p.x - p.y, p.y + p.z) - 2.0 * s) <= 0.0 )
    {
        d = max(d_hex, d_tri);
        g = (d_tri > d_hex) ? vec3(0.57735027) : vec3(0.57735027, 0.57735027, -0.57735027);
    }
    else
    {
        float t    = max(c - 0.5 * min(p.x - p.y, p.y + p.z), 0.0);
        vec3  v    = (p.x - p.y <= p.y + p.z)
                   ? vec3(p.x - s + t, p.y - c - t, p.z - c)
                   : vec3(p.x - s, p.y - c + t, p.z - c + t);
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
