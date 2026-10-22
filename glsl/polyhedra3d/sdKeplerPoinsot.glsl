// ============================================================================
// sdSmallStellatedDodecahedron, sdGreatDodecahedron, sdGreatStellatedDodecahedron
// & sdg... — Exact 0/1-Sqrt 3D Dodecahedral Kepler-Poinsot Regular Star Polyhedra
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in face prisms, 1 sqrt in ridges/valleys)
// @novelty Exact 0/1-sqrt Euclidean SDFs and analytical gradients for the three dodecahedral
//          Kepler-Poinsot regular star polyhedra ({5/2, 5}, {5, 5/2}, {5/2, 3}) via the I_h
//          Order-Statistic Total-Ordering Theorem (ax + bz >= -ax + bz >= by + az >= -by + az
//          for all p in C_{I_h}) and exact coordinate swizzles -n3.yzx and -n3.zxy for the
//          oblique ridge/valley unit tangents.
//
// Parameter r: perpendicular distance from the origin to the star face planes.
// ============================================================================

#ifndef SDF_GEMS_SD_KEPLER_POINSOT_GLSL
#define SDF_GEMS_SD_KEPLER_POINSOT_GLSL

// 1. Small Stellated Dodecahedron {5/2, 5} (1st stellation of the dodecahedron, k = 2)
float sdSmallStellatedDodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    const vec2 n = vec2(0.52573111, 0.85065081);
    float df  = -n.x * p.x + n.y * p.z - r;
    float dn3 = dot(p, n3);
    if( n.x * df >= -p.x && n.x * df <= dn3 )
        return df;

    float y3 = 0.44902798 * r, z3 = 1.17557050 * r;
    vec3  q3 = p - vec3(0.0, y3, z3);
    vec3  e  = -n3.yzx; // (0.80901699, -0.30901699, 0.5)
    vec3  v  = (df < 0.0)
             ? vec3(p.x, max(q3.y, 0.0), q3.z)
             : q3 - clamp(dot(q3, e), 0.0, 1.45308506 * r) * e;
    return ((df < 0.0) ? -1.0 : 1.0) * length(v);
}

// 2. Great Dodecahedron {5, 5/2} (2nd stellation of the dodecahedron, k = 3)
float sdGreatDodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    const vec2 n = vec2(0.52573111, 0.85065081);
    float df  = n.y * p.y + n.x * p.z - r;
    float dn3 = dot(p, n3);
    if( n.x * df >= -dn3 && n.y * df <= p.y )
        return df;

    float y3 = 0.44902798 * r, z3 = 1.17557050 * r, z5 = 1.90211303 * r;
    vec3  q3 = p - vec3(0.0, y3, z3);
    vec3  e  = -n3.yzx; // (0.80901699, -0.30901699, 0.5)
    vec3  v  = (df < 0.0)
             ? q3 - clamp(dot(q3, e), 0.0, 1.45308506 * r) * e
             : vec3(max(p.x - z3, 0.0), p.y, p.z - z5);
    return ((df < 0.0) ? -1.0 : 1.0) * length(v);
}

// 3. Great Stellated Dodecahedron {5/2, 3} (3rd stellation of the dodecahedron, k = 4)
float sdGreatStellatedDodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    const vec2 n = vec2(0.52573111, 0.85065081);
    float df  = -n.y * p.y + n.x * p.z - r;
    float dn3 = dot(p, n3);
    if( n.y * df >= -p.y && n.y * df <= dn3 )
        return df;

    float x5 = 1.17557050 * r, z5 = 1.90211303 * r;
    vec3  q5 = p - vec3(x5, 0.0, z5);
    vec3  e  = -n3.zxy; // (-0.30901699, 0.5, 0.80901699)
    vec3  v  = (df < 0.0)
             ? vec3(max(q5.x, 0.0), p.y, q5.z)
             : q5 - clamp(dot(q5, e), 0.0, 3.80422607 * r) * e;
    return ((df < 0.0) ? -1.0 : 1.0) * length(v);
}

// 4. Analytical Gradients sdgSmallStellatedDodecahedron, sdgGreatDodecahedron, sdgGreatStellatedDodecahedron
vec4 sdgSmallStellatedDodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    vec3 s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); bool r0 = (dot(p, n3) < 0.0); if( r0 ) p -= 2.0 * dot(p, n3) * n3;
    vec3 s1 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); bool r1 = (dot(p, n3) < 0.0); if( r1 ) p -= 2.0 * dot(p, n3) * n3;
    vec3 s2 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); bool r2 = (dot(p, n3) < 0.0); if( r2 ) p -= 2.0 * dot(p, n3) * n3;
    vec2 s3 = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    p.xy = abs(p.xy);

    const vec2 n = vec2(0.52573111, 0.85065081);
    float df  = -n.x * p.x + n.y * p.z - r;
    float dn3 = dot(p, n3);
    float d;
    vec3  g;
    if( n.x * df >= -p.x && n.x * df <= dn3 )
    {
        d = df;
        g = vec3(-n.x, 0.0, n.y);
    }
    else
    {
        float y3   = 0.44902798 * r, z3 = 1.17557050 * r;
        vec3  q3   = p - vec3(0.0, y3, z3);
        vec3  e    = -n3.yzx;
        vec3  v    = (df < 0.0)
                   ? vec3(p.x, max(q3.y, 0.0), q3.z)
                   : q3 - clamp(dot(q3, e), 0.0, 1.45308506 * r) * e;
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        float sm   = (df < 0.0) ? -1.0 : 1.0;
        d = sm * v2 * invV;
        g = (sm * invV) * v;
    }

    g.xy *= s3;
    if( r2 ) g -= 2.0 * dot(g, n3) * n3; g *= s2;
    if( r1 ) g -= 2.0 * dot(g, n3) * n3; g *= s1;
    if( r0 ) g -= 2.0 * dot(g, n3) * n3; g *= s0;
    return vec4(d, g);
}

vec4 sdgGreatDodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    vec3 s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); bool r0 = (dot(p, n3) < 0.0); if( r0 ) p -= 2.0 * dot(p, n3) * n3;
    vec3 s1 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); bool r1 = (dot(p, n3) < 0.0); if( r1 ) p -= 2.0 * dot(p, n3) * n3;
    vec3 s2 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); bool r2 = (dot(p, n3) < 0.0); if( r2 ) p -= 2.0 * dot(p, n3) * n3;
    vec2 s3 = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    p.xy = abs(p.xy);

    const vec2 n = vec2(0.52573111, 0.85065081);
    float df  = n.y * p.y + n.x * p.z - r;
    float dn3 = dot(p, n3);
    float d;
    vec3  g;
    if( n.x * df >= -dn3 && n.y * df <= p.y )
    {
        d = df;
        g = vec3(0.0, n.y, n.x);
    }
    else
    {
        float y3   = 0.44902798 * r, z3 = 1.17557050 * r, z5 = 1.90211303 * r;
        vec3  q3   = p - vec3(0.0, y3, z3);
        vec3  e    = -n3.yzx;
        vec3  v    = (df < 0.0)
                   ? q3 - clamp(dot(q3, e), 0.0, 1.45308506 * r) * e
                   : vec3(max(p.x - z3, 0.0), p.y, p.z - z5);
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        float sm   = (df < 0.0) ? -1.0 : 1.0;
        d = sm * v2 * invV;
        g = (sm * invV) * v;
    }

    g.xy *= s3;
    if( r2 ) g -= 2.0 * dot(g, n3) * n3; g *= s2;
    if( r1 ) g -= 2.0 * dot(g, n3) * n3; g *= s1;
    if( r0 ) g -= 2.0 * dot(g, n3) * n3; g *= s0;
    return vec4(d, g);
}

vec4 sdgGreatStellatedDodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    vec3 s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); bool r0 = (dot(p, n3) < 0.0); if( r0 ) p -= 2.0 * dot(p, n3) * n3;
    vec3 s1 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); bool r1 = (dot(p, n3) < 0.0); if( r1 ) p -= 2.0 * dot(p, n3) * n3;
    vec3 s2 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); bool r2 = (dot(p, n3) < 0.0); if( r2 ) p -= 2.0 * dot(p, n3) * n3;
    vec2 s3 = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    p.xy = abs(p.xy);

    const vec2 n = vec2(0.52573111, 0.85065081);
    float df  = -n.y * p.y + n.x * p.z - r;
    float dn3 = dot(p, n3);
    float d;
    vec3  g;
    if( n.y * df >= -p.y && n.y * df <= dn3 )
    {
        d = df;
        g = vec3(0.0, -n.y, n.x);
    }
    else
    {
        float x5   = 1.17557050 * r, z5 = 1.90211303 * r;
        vec3  q5   = p - vec3(x5, 0.0, z5);
        vec3  e    = -n3.zxy;
        vec3  v    = (df < 0.0)
                   ? vec3(max(q5.x, 0.0), p.y, q5.z)
                   : q5 - clamp(dot(q5, e), 0.0, 3.80422607 * r) * e;
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        float sm   = (df < 0.0) ? -1.0 : 1.0;
        d = sm * v2 * invV;
        g = (sm * invV) * v;
    }

    g.xy *= s3;
    if( r2 ) g -= 2.0 * dot(g, n3) * n3; g *= s2;
    if( r1 ) g -= 2.0 * dot(g, n3) * n3; g *= s1;
    if( r0 ) g -= 2.0 * dot(g, n3) * n3; g *= s0;
    return vec4(d, g);
}

#endif // SDF_GEMS_SD_KEPLER_POINSOT_GLSL
