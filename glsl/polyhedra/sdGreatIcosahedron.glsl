// ============================================================================
// sdGreatIcosahedron & sdgGreatIcosahedron — Exact 0/1-Sqrt 3D Great Icosahedron {3, 5/2}
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in face prisms, 1 sqrt in ridges/valleys)
// @novelty Exact 0/1-sqrt Euclidean SDF and analytical gradient for the 4th Kepler-Poinsot
//          regular star polyhedron, the Great Icosahedron {3, 5/2} (7th stellation of the
//          icosahedron with 180 visible outer triangular facets), via I_h Weyl chamber folding,
//          collapsing 180 facets into the union of two planes d_f = min(d_A, d_B) <= 0 with
//          N_A . N_B = -1/3 and algebraic saddle valley tangent (-(phi-1), sqrt(5), phi)/(2*sqrt(2)).
//
// Parameter r: perpendicular distance from the origin to the 20 equilateral star face planes.
// ============================================================================

#ifndef SDF_GEMS_SD_GREAT_ICOSAHEDRON_GLSL
#define SDF_GEMS_SD_GREAT_ICOSAHEDRON_GLSL

float sdGreatIcosahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    const float inv3 = 0.57735027;
    const vec2  nB   = vec2(0.93417236, 0.35682209);
    float dA  = (-p.x - p.y + p.z) * inv3 - r;
    float dB  = nB.x * p.x + nB.y * p.z - r;
    float df  = min(dA, dB);
    float dn3 = dot(p, n3);

    bool inFA_ext = (df > 0.0 && dn3 >= nB.x * dA);
    bool inFB_ext = (df > 0.0 && p.x >= nB.x * dB);
    if( (dA <= 0.0 && p.y + inv3 * dA >= 0.0 && dA + 3.0 * dB >= 0.0) ||
        (dB <= 0.0 && 3.0 * dA + dB >= 0.0) ||
        ((dA <= dB) ? inFA_ext : inFB_ext) )
        return df;

    if( df <= 0.0 )
    {
        vec3 qm = p - vec3(0.29586960 * r, 0.0, 2.02792041 * r);
        const vec3 ev1 = vec3(0.70710678, 0.0, 0.70710678);
        const vec3 ev2 = vec3(-0.21850801, 0.79056942, 0.57206140);
        vec3 v1 = qm - clamp(dot(qm, ev1), 0.0, 3.54493482 * r) * ev1;
        vec3 v2 = qm - clamp(dot(qm, ev2), 0.0, 1.35404463 * r) * ev2;
        return -sqrt(min(dot(v1, v1), dot(v2, v2)));
    }

    float y3 = 1.07046627 * r, z3 = 2.80251708 * r;
    vec3  q3 = p - vec3(0.0, y3, z3);
    vec3  eA = -n3.yzx;
    vec3  vA = q3 - clamp(dot(q3, eA), 0.0, 3.46410162 * r) * eA;
    vec3  vB = vec3(p.x, max(q3.y, 0.0), q3.z);
    float d2A = inFA_ext ? (dA * dA) : dot(vA, vA);
    float d2B = inFB_ext ? (dB * dB) : dot(vB, vB);
    return sqrt(min(d2A, d2B));
}

vec4 sdgGreatIcosahedron( vec3 p, float r )
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

    const float inv3 = 0.57735027;
    const vec3  nA   = vec3(-0.57735027, -0.57735027, 0.57735027);
    const vec3  nB   = vec3(0.93417236, 0.0, 0.35682209);
    float dA  = dot(p, nA) - r;
    float dB  = dot(p, nB) - r;
    float df  = min(dA, dB);
    float dn3 = dot(p, n3);

    bool inFA_ext = (df > 0.0 && dn3 >= nB.x * dA);
    bool inFB_ext = (df > 0.0 && p.x >= nB.x * dB);
    float d;
    vec3  g;
    if( (dA <= 0.0 && p.y + inv3 * dA >= 0.0 && dA + 3.0 * dB >= 0.0) ||
        (dB <= 0.0 && 3.0 * dA + dB >= 0.0) ||
        ((dA <= dB) ? inFA_ext : inFB_ext) )
    {
        d = df;
        g = (dA <= dB) ? nA : nB;
    }
    else
    {
        float y3 = 1.07046627 * r, z3 = 2.80251708 * r;
        vec3  q3 = p - vec3(0.0, y3, z3);
        vec3  eA = -n3.yzx;
        vec3  vA = inFA_ext ? (dA * nA) : (q3 - clamp(dot(q3, eA), 0.0, 3.46410162 * r) * eA);
        vec3  vB = inFB_ext ? (dB * nB) : vec3(p.x, max(q3.y, 0.0), q3.z);

        vec3  qm = p - vec3(0.29586960 * r, 0.0, 2.02792041 * r);
        const vec3 ev1 = vec3(0.70710678, 0.0, 0.70710678);
        const vec3 ev2 = vec3(-0.21850801, 0.79056942, 0.57206140);
        vec3  vv1 = qm - clamp(dot(qm, ev1), 0.0, 3.54493482 * r) * ev1;
        vec3  vv2 = qm - clamp(dot(qm, ev2), 0.0, 1.35404463 * r) * ev2;

        vec3  v1   = (df < 0.0) ? vv1 : vA;
        vec3  v2   = (df < 0.0) ? vv2 : vB;
        float d1   = dot(v1, v1), d2 = dot(v2, v2);
        vec3  v    = (d1 <= d2) ? v1 : v2;
        float v2m  = min(d1, d2);
        float invV = inversesqrt(max(v2m, 1e-30));
        float sm   = (df < 0.0) ? -1.0 : 1.0;
        d = sm * v2m * invV;
        g = (sm * invV) * v;
    }

    g.xy *= s3;
    if( r2 ) g -= 2.0 * dot(g, n3) * n3; g *= s2;
    if( r1 ) g -= 2.0 * dot(g, n3) * n3; g *= s1;
    if( r0 ) g -= 2.0 * dot(g, n3) * n3; g *= s0;
    return vec4(d, g);
}

#endif // SDF_GEMS_SD_GREAT_ICOSAHEDRON_GLSL
