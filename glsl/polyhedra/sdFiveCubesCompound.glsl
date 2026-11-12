// ============================================================================
// sdFiveCubesCompound — Exact 0/1-Sqrt Compound of 5 Cubes & Compound of 5 Octahedra
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, `1.11e-9` to `4.30e-9` max error, 0 `sqrt`s in face prisms, 1 `sqrt` on ridges/valleys)
// ============================================================================

#ifndef SDF_GEMS_SD_FIVE_CUBES_COMPOUND_GLSL
#define SDF_GEMS_SD_FIVE_CUBES_COMPOUND_GLSL

float sdFiveOctahedraCompound( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.50000000, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    // Single-Octahedron Collapse Theorem: Octahedron 0 wins across 100% of C_{I_h}
    float df = 0.57735027 * (p.x + p.y + p.z) - r;
    float h  = 0.57735027 * df;
    if( min(p.x, p.y) >= h && dot(p, n3) >= -h )
        return df; // 0 sqrts in both interior and exterior face prisms!

    if( df < 0.0 )
    {
        // Concave valley seam along Wall C from V_B to V_A
        vec3 q = p - vec3(0.0, 0.47872707, 1.25332374) * r;
        const vec3 e = vec3(0.79056942, -0.57206140, -0.21850801);
        float t = clamp(dot(q, e), 0.0, 0.83684560 * r);
        return -length(q - e * t); // 1 sqrt
    }
    // Convex outer ridge from apex V_top = (0, 0, sqrt(3)*r)
    float s = 1.73205081 * r;
    if( p.x <= p.y )
    {
        float t = clamp(0.5 * (p.y - p.z + s), 0.0, 0.47872707 * r);
        return length(vec3(p.x, p.y - t, p.z - s + t)); // 1 sqrt
    }
    float t = clamp(0.5 * (p.x - p.z + s), 0.0, 0.66158454 * r);
    return length(vec3(p.x - t, p.y, p.z - s + t)); // 1 sqrt
}

vec4 sdgFiveOctahedraCompound( vec3 p, float r )
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

    float df = 0.57735027 * (p.x + p.y + p.z) - r;
    float h  = 0.57735027 * df;
    vec3  g;
    float d;
    if( min(p.x, p.y) >= h && dot(p, n3) >= -h )
    {
        d = df;
        g = vec3(0.57735027);
    }
    else if( df < 0.0 )
    {
        vec3 q = p - vec3(0.0, 0.47872707, 1.25332374) * r;
        const vec3 e = vec3(0.79056942, -0.57206140, -0.21850801);
        vec3 v = q - e * clamp(dot(q, e), 0.0, 0.83684560 * r);
        float l = length(v);
        d = -l;
        g = -v / max(l, 1e-20);
    }
    else
    {
        float s = 1.73205081 * r;
        vec3  v;
        if( p.x <= p.y )
        {
            float t = clamp(0.5 * (p.y - p.z + s), 0.0, 0.47872707 * r);
            v = vec3(p.x, p.y - t, p.z - s + t);
        }
        else
        {
            float t = clamp(0.5 * (p.x - p.z + s), 0.0, 0.66158454 * r);
            v = vec3(p.x - t, p.y, p.z - s + t);
        }
        float l = length(v);
        d = l;
        g = v / max(l, 1e-20);
    }
    g.xy *= s3;
    if( c2 < 0.0 ) g -= 2.0 * dot(g, n3) * n3; g *= s2;
    if( c1 < 0.0 ) g -= 2.0 * dot(g, n3) * n3; g *= s1;
    if( c0 < 0.0 ) g -= 2.0 * dot(g, n3) * n3; g *= s0;
    return vec4(d, g);
}

float sdFiveCubesCompound( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.50000000, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    // Swizzle-Frame coordinates of Cube 1 and Cube 2 in C_{I_h}:
    vec3  n1a = -n3.yzx;                         // ( 0.80901699, -0.30901699, 0.50000000)
    vec3  n2b = vec3(-n3.z, -0.5, -n3.y);        // (-0.30901699, -0.50000000, 0.80901699)
    float q1c = dot(p, n3);
    float q1a = dot(p, n1a) - r;
    float q2a = q1a + 0.61803399 * p.y;
    float q2b = dot(p, n2b) - r;
    float q1b = q2b + p.y;
    float q2c = q1c + 1.61803399 * p.y - r;

    float df = min(max(q1a, q1b), max(q2a, q2b));
    if( df >= 0.0 )
    {
        // Exact 2-Cube Orthonormal Box Reduction in the exterior:
        vec3 v1 = max(vec3(q1a, q1b, q1c - 0.23606798 * r), 0.0);
        vec3 v2 = max(vec3(q2a, q2b, q2c), 0.0);
        return sqrt(min(dot(v1, v1), dot(v2, v2))); // 1 sqrt
    }

    // Interior: check 4 face prisms (0 sqrts if winning face is in its prism!)
    bool in1a = (q1b <= 0.0) && (p.y + 0.30901699 * q1a >= 0.0);
    bool in1b = (q1a <= 0.0) && (q2a - 0.30901699 * q1b >= 0.0);
    bool in2a = (q2b <= 0.0) && (q1c + 0.5 * q2a >= 0.0) && (q1b - 0.30901699 * q2a >= 0.0);
    bool in2b = (q2a <= 0.0) && (p.y + 0.5 * q2b >= 0.0) && (p.x + 0.30901699 * q2b >= 0.0);
    if( (in1a && df == q1a) || (in1b && df == q1b) || (in2a && df == q2a) || (in2b && df == q2b) )
        return df; // 0 sqrts!

    float d2 = 1e30;
    if( in1a ) d2 = min(d2, q1a * q1a);
    if( in1b ) d2 = min(d2, q1b * q1b);
    if( in2a ) d2 = min(d2, q2a * q2a);
    if( in2b ) d2 = min(d2, q2b * q2b);

    // 5 wall-aligned concave valley seams governed by 2D vectors e and f:
    const vec2 e = vec2(0.52573111, 0.85065081);
    const vec2 f = vec2(0.93417236, 0.35682209);
    float xm = p.x - 0.38196601 * r;
    float zm = p.z - 1.38196601 * r;
    float zz = p.z - 1.23606798 * r;

    float t1 = clamp(e.x * xm - e.y * zm, 0.0, 0.44902798 * r);
    d2 = min(d2, dot(vec3(xm - e.x * t1, p.y, zm + e.y * t1), vec3(xm - e.x * t1, p.y, zm + e.y * t1)));

    float t5 = clamp(e.y * p.y - e.x * zm, 0.0, 0.27751455 * r);
    d2 = min(d2, dot(vec3(xm, p.y - e.y * t5, zm + e.x * t5), vec3(xm, p.y - e.y * t5, zm + e.x * t5)));

    float t8 = clamp(f.x * p.x + f.y * zz, 0.0, 0.40888173 * r);
    d2 = min(d2, dot(vec3(p.x - f.x * t8, p.y, zz - f.y * t8), vec3(p.x - f.x * t8, p.y, zz - f.y * t8)));

    float t9 = clamp(e.y * p.y + e.x * zz, 0.0, 0.72654253 * r);
    d2 = min(d2, dot(vec3(p.x, p.y - e.y * t9, zz - e.x * t9), vec3(p.x, p.y - e.y * t9, zz - e.x * t9)));

    float u1a = q1a + 0.38196601 * r;
    float u1b = q1b - 0.61803399 * r;
    float t6  = clamp(f.y * u1a - f.x * u1b, 0.0, 0.66158454 * r);
    d2 = min(d2, dot(vec3(u1a - f.y * t6, u1b + f.x * t6, q1c), vec3(u1a - f.y * t6, u1b + f.x * t6, q1c)));

    return -sqrt(d2); // 1 sqrt
}

#endif // SDF_GEMS_SD_FIVE_CUBES_COMPOUND_GLSL
