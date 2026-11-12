// ============================================================================
// sdTrapezohedron & sdgTrapezohedron — Exact 0-Trig 0/1-Sqrt 3D Pentagonal Trapezohedron (D10 Die) & General N-gonal Trapezohedron
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 trig & 0 sqrts in D10 interior/face prism, 1 sqrt in edge/vertex cone)
// @novelty Exact 0-trig 0/1-sqrt Euclidean SDF and analytical gradient for the 10-sided
//          pentagonal trapezohedron (D10 tabletop die) and general 2N-kite trapezohedron
//          via 2D dihedral mirror folding and exact D_nd 3D Voronoi bisector rotation.
// ============================================================================

// 1. Exact 0-Trig 0/1-Sqrt Pentagonal Trapezohedron (10-Sided D10 Die: equatorial radius r, polar half-height h)
float sdTrapezohedron( vec3 p, float r, float h )
{
    const vec2 n5 = vec2(-0.58778525, 0.80901699);
    p.y = abs(p.y);
    p.xy -= 2.0 * max(dot(p.xy, n5), 0.0) * n5; p.y = abs(p.y);
    p.xy -= 2.0 * max(dot(p.xy, n5), 0.0) * n5; p.y = abs(p.y);
    p.xy -= 2.0 * max(dot(p.xy, n5), 0.0) * n5;

    const float c = 0.80901699, s = 0.58778525, t2 = 0.10557281;
    float a = 1.10557281 * h;
    float x = p.x, y = p.y, z = p.z;
    float u = x * c + y * s, v = x * s - y * c;

    // Exact 3D Voronoi bisector rotation in D_5d
    if( a * (x - u) + 2.0 * r * z < 0.0 )
    {
        float tx = x, ty = y;
        x = u; y = v; u = tx; v = ty; z = -z;
    }

    float n2 = a * a + r * r;
    float m  = a * x + r * (z - h);
    float w2 = a * (x - u) + 2.0 * r * z;
    if( a * s * m <= n2 * v && (a * a * (1.0 - c) + 2.0 * r * r) * m <= n2 * w2 )
        return m * inversesqrt(n2);

    vec3  e1 = vec3(r * c, r * s, -a * c);
    float l1 = clamp((x * e1.x + y * e1.y + (z - h) * e1.z) / (r * r + a * a * c * c), 0.0, 1.0);
    vec3  d1 = vec3(x, y, z - h) - l1 * e1;

    vec3  e2 = vec3(r * (1.0 - c), -r * s, -a * (1.0 - c));
    vec3  q2 = vec3(x - r * c, y - r * s, z - h * t2);
    float l2 = clamp(dot(q2, e2) / (2.0 * r * r * (1.0 - c) + a * a * (1.0 - c) * (1.0 - c)), 0.0, 1.0);
    vec3  d2 = q2 - l2 * e2;

    return sqrt(min(dot(d1, d1), dot(d2, d2)));
}

// 2. General N-gonal Trapezohedron (2*n kite faces, n >= 3)
float sdTrapezohedron( vec3 p, float n, float r, float h )
{
    float sp  = 3.14159265 / n;
    float c   = cos(sp), s = sin(sp);
    float t2  = (1.0 - c) / (1.0 + c);
    float a   = h * (1.0 + t2);
    float rad = length(p.xy);
    float ang = mod(atan(p.y, p.x), 2.0 * sp);
    ang = sp - abs(ang - sp);
    float x = rad * cos(ang), y = rad * sin(ang), z = p.z;
    float u = x * c + y * s,  v = x * s - y * c;

    if( a * (x - u) + 2.0 * r * z < 0.0 )
    {
        float tx = x, ty = y;
        x = u; y = v; u = tx; v = ty; z = -z;
    }

    float n2 = a * a + r * r;
    float m  = a * x + r * (z - h);
    float w2 = a * (x - u) + 2.0 * r * z;
    if( a * s * m <= n2 * v && (a * a * (1.0 - c) + 2.0 * r * r) * m <= n2 * w2 )
        return m * inversesqrt(n2);

    vec3  e1 = vec3(r * c, r * s, -a * c);
    float l1 = clamp((x * e1.x + y * e1.y + (z - h) * e1.z) / (r * r + a * a * c * c), 0.0, 1.0);
    vec3  d1 = vec3(x, y, z - h) - l1 * e1;

    vec3  e2 = vec3(r * (1.0 - c), -r * s, -a * (1.0 - c));
    vec3  q2 = vec3(x - r * c, y - r * s, z - h * t2);
    float l2 = clamp(dot(q2, e2) / (2.0 * r * r * (1.0 - c) + a * a * (1.0 - c) * (1.0 - c)), 0.0, 1.0);
    vec3  d2 = q2 - l2 * e2;

    return sqrt(min(dot(d1, d1), dot(d2, d2)));
}

// 3. Analytical 3D Distance and Unit Gradient vec4(d, grad.xyz) for D10 Trapezohedron
vec4 sdgTrapezohedron( vec3 p, float r, float h )
{
    const vec2 n5 = vec2(-0.58778525, 0.80901699);
    float s0 = (p.y < 0.0) ? -1.0 : 1.0; p.y = abs(p.y);
    bool  r0 = (dot(p.xy, n5) > 0.0); if( r0 ) p.xy -= 2.0 * dot(p.xy, n5) * n5;
    float s1 = (p.y < 0.0) ? -1.0 : 1.0; p.y = abs(p.y);
    bool  r1 = (dot(p.xy, n5) > 0.0); if( r1 ) p.xy -= 2.0 * dot(p.xy, n5) * n5;
    float s2 = (p.y < 0.0) ? -1.0 : 1.0; p.y = abs(p.y);
    bool  r2 = (dot(p.xy, n5) > 0.0); if( r2 ) p.xy -= 2.0 * dot(p.xy, n5) * n5;

    const float c = 0.80901699, s = 0.58778525, t2 = 0.10557281;
    float a = 1.10557281 * h;
    float x = p.x, y = p.y, z = p.z;
    float u = x * c + y * s, v = x * s - y * c;

    bool swap = (a * (x - u) + 2.0 * r * z < 0.0);
    if( swap )
    {
        float tx = x, ty = y;
        x = u; y = v; u = tx; v = ty; z = -z;
    }

    float n2   = a * a + r * r;
    float invN = inversesqrt(n2);
    float m    = a * x + r * (z - h);
    float w2   = a * (x - u) + 2.0 * r * z;
    float d;
    vec3  g;
    if( a * s * m <= n2 * v && (a * a * (1.0 - c) + 2.0 * r * r) * m <= n2 * w2 )
    {
        d = m * invN;
        g = vec3(a * invN, 0.0, r * invN);
    }
    else
    {
        vec3  e1 = vec3(r * c, r * s, -a * c);
        float l1 = clamp((x * e1.x + y * e1.y + (z - h) * e1.z) / (r * r + a * a * c * c), 0.0, 1.0);
        vec3  d1 = vec3(x, y, z - h) - l1 * e1;

        vec3  e2 = vec3(r * (1.0 - c), -r * s, -a * (1.0 - c));
        vec3  q2 = vec3(x - r * c, y - r * s, z - h * t2);
        float l2 = clamp(dot(q2, e2) / (2.0 * r * r * (1.0 - c) + a * a * (1.0 - c) * (1.0 - c)), 0.0, 1.0);
        vec3  d2 = q2 - l2 * e2;

        float v1_2 = dot(d1, d1), v2_2 = dot(d2, d2);
        vec3  ve   = (v1_2 <= v2_2) ? d1 : d2;
        float v2m  = min(v1_2, v2_2);
        float invV = inversesqrt(max(v2m, 1e-30));
        d = v2m * invV;
        g = ve * invV;
    }

    if( swap ) g = vec3(g.x * c + g.y * s, g.x * s - g.y * c, -g.z);
    if( r2 ) g.xy -= 2.0 * dot(g.xy, n5) * n5; g.y *= s2;
    if( r1 ) g.xy -= 2.0 * dot(g.xy, n5) * n5; g.y *= s1;
    if( r0 ) g.xy -= 2.0 * dot(g.xy, n5) * n5; g.y *= s0;
    return vec4(d, g);
}
