// ============================================================================
// sdgOctahedron — 1-Sqrt Polyhedral Analytical 3D Gradients via Weyl Decomposition
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1`, `1` sqrt, branchless Weyl-chamber gradient formula)
// ============================================================================

// Helper: sort q = abs(p) into x >= y >= z >= 0 and invert the permutation on a Weyl-chamber vector g.
void _sdgWeylSort( vec3 p, out vec3 q, out vec3 s0, out bvec3 sw )
{
    s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    q  = abs(p);
    sw.x = (q.x < q.y); if( sw.x ) q.xy = q.yx;
    sw.y = (q.y < q.z); if( sw.y ) q.yz = q.zy;
    sw.z = (q.x < q.y); if( sw.z ) q.xy = q.yx;
}

vec3 _sdgWeylUnsort( vec3 g, vec3 s0, bvec3 sw )
{
    if( sw.z ) g.xy = g.yx;
    if( sw.y ) g.yz = g.zy;
    if( sw.x ) g.xy = g.yx;
    return s0 * g;
}

// 1. Regular Octahedron (|x| + |y| + |z| <= s)
//    In sorted Weyl chamber x >= y >= z >= 0:
//      m = x + y + z - s,  a = max(m - 3z, 0),  b = max(x - y - s, 0),
//      d = sign(m) * sqrt((2*m^2 + a^2 + 3*b^2) / 6),
//      g_sorted = (2*m + a + 3*b,  2*m + a - 3*b,  2*m - 2*a) / (6*d).
vec4 sdgOctahedron( vec3 p, float s )
{
    vec3  q, s0;
    bvec3 sw;
    _sdgWeylSort(p, q, s0, sw);

    float m  = q.x + q.y + q.z - s;
    float a  = max(m - 3.0 * q.z, 0.0);
    float b  = max(q.x - q.y - s, 0.0);
    float sm = (m < 0.0) ? -1.0 : 1.0;
    float d  = sm * sqrt((2.0 * m * m + a * a + 3.0 * b * b) * (1.0 / 6.0));
    vec3  g  = (abs(d) > 0.0)
             ? (vec3(2.0 * m + a + 3.0 * b, 2.0 * m + a - 3.0 * b, 2.0 * m - 2.0 * a) / (6.0 * d))
             : vec3(0.57735027);
    return vec4(d, _sdgWeylUnsort(g, s0, sw));
}

// 2. Catalan Rhombic Dodecahedron (max(|x|+|y|, |y|+|z|, |z|+|x|) <= s)
//    In sorted Weyl chamber x >= y >= z >= 0:
//      m = x + y - s,  a = max(x - y + 2z - s, 0),
//      b0 = y + z - x + 0.25*s,  b = sign(b0) * max(abs(b0) - 0.75*s, 0),
//      d = sign(m) * sqrt((3*m^2 + a^2 + 2*b^2) / 6),
//      g_sorted = (3*m + a - 2*b,  3*m - a + 2*b,  2*a + 2*b) / (6*d).
vec4 sdgRhombicDodecahedron( vec3 p, float s )
{
    vec3  q, s0;
    bvec3 sw;
    _sdgWeylSort(p, q, s0, sw);

    float m  = q.x + q.y - s;
    float a  = max(q.x - q.y + 2.0 * q.z - s, 0.0);
    float b0 = q.y + q.z - q.x + 0.25 * s;
    float b  = ((b0 < 0.0) ? -1.0 : 1.0) * max(abs(b0) - 0.75 * s, 0.0);
    float sm = (m < 0.0) ? -1.0 : 1.0;
    float d  = sm * sqrt((3.0 * m * m + a * a + 2.0 * b * b) * (1.0 / 6.0));
    vec3  g  = (abs(d) > 0.0)
             ? (vec3(3.0 * m + a - 2.0 * b, 3.0 * m - a + 2.0 * b, 2.0 * a + 2.0 * b) / (6.0 * d))
             : vec3(0.70710678, 0.70710678, 0.0);
    return vec4(d, _sdgWeylUnsort(g, s0, sw));
}
