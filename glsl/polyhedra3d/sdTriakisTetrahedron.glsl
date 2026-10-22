// ============================================================================
// sdTriakisTetrahedron & sdgTriakisTetrahedron — Exact 0/1-Sqrt 3D Catalan Triakis Tetrahedron
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & face prism, 1 sqrt in edge/vertex cone)
// ============================================================================

// 1. Exact 0/1-Sqrt Catalan Triakis Tetrahedron (12 isosceles triangular faces, dual of Truncated Tetrahedron)
float sdTriakisTetrahedron( vec3 p, float s )
{
    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;
    p.z *= sig;

    float m = 3.0 * p.x + p.y + p.z - 3.0 * s;
    if( m <= 5.5 * min(p.x - p.y, p.y + p.z) )
        return m * 0.30151134; // 1.0 / sqrt(11.0): 0 sqrts in interior & face prism!

    float t1 = clamp((6.0 * s - p.x - p.y + 4.0 * p.z) * (1.0 / 18.0), 0.0, 0.4 * s);
    float u2 = min(0.5 * (p.y - p.z), s);
    vec3  d  = (p.x - p.y <= p.y + p.z)
             ? vec3(p.x - s + t1, p.y - s + t1, p.z + s - 4.0 * t1)
             : vec3(p.x - s, p.y - u2, p.z + u2);
    return length(d);
}

// 2. Ultra-Compact Branchless Pythagorean Weyl-Basis Form (7 lines, 0 sorting branches, 1 sqrt)
float sdTriakisTetrahedronBranchless( vec3 p, float s )
{
    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
    p = abs(p);
    float x = max(p.x, max(p.y, p.z)), z = min(p.x, min(p.y, p.z));
    float y = p.x + p.y + p.z - x - z; z *= sig;
    float m = 3.0 * x + y + z - 3.0 * s, a = max(2.0 * m - 11.0 * min(x - y, y + z), 0.0);
    float b = (x - y <= y + z) ? max(abs(4.0 * z - x - y + 2.4 * s) - 3.6 * s, 0.0)
                               : max(3.0 * (y - z) - 6.0 * s, 0.0);
    return ((m < 0.0) ? -1.0 : 1.0) * sqrt((18.0 * m * m + a * a + 11.0 * b * b) / 198.0);
}

// 3. Analytical 3D Distance and Unit Gradient vec4(d, grad.xyz) (0/1-inversesqrt)
vec4 sdgTriakisTetrahedron( vec3 p, float s )
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

    float m = 3.0 * p.x + p.y + p.z - 3.0 * s;
    float d;
    vec3  g;
    if( m <= 5.5 * min(p.x - p.y, p.y + p.z) )
    {
        d = m * 0.30151134;
        g = vec3(0.90453403, 0.30151134, 0.30151134);
    }
    else
    {
        float t1   = clamp((6.0 * s - p.x - p.y + 4.0 * p.z) * (1.0 / 18.0), 0.0, 0.4 * s);
        float u2   = min(0.5 * (p.y - p.z), s);
        vec3  v    = (p.x - p.y <= p.y + p.z)
                   ? vec3(p.x - s + t1, p.y - s + t1, p.z + s - 4.0 * t1)
                   : vec3(p.x - s, p.y - u2, p.z + u2);
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
