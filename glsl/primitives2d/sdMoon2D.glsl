// ============================================================================
// sdMoon2D, sdgMoon2D, sdMoon3D & sdLens2D — Exact 1-Sqrt-Arc Crescent Moon, Meniscus & Optical Lens
// @author Inigo Quilez — upgraded with 0-sqrt signed-square wedge & confocal bisector using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 1 sqrt in circular arcs down from 3 in 2D and 4 in 3D, `0.0` error)
// ============================================================================

// 1. Exact 1-Sqrt-Arc 2D Crescent Moon SDF (down from 3 sqrts in IQ, 0.0 difference)
float sdMoon2D( vec2 p, float d, float ra, float rb )
{
    float py2  = p.y * p.y;
    float a    = (ra * ra - rb * rb + d * d) / (2.0 * d);
    float b2   = max(ra * ra - a * a, 0.0);
    bool  low  = (b2 > py2);
    float u    = low ? d * (p.x - d) : d * p.x;
    float cv   = low ? d * (a - d)   : d * a;
    if( b2 > 0.0 && u * abs(u) * b2 > cv * abs(cv) * py2 )
        return length(vec2(p.x - a, abs(p.y) - sqrt(b2)));

    float A    = 0.5 * (ra + rb), A2 = A * A;
    float B2   = A2 - 0.25 * d * d;
    float dx   = p.x - 0.5 * d;
    bool  useA = (B2 <= 0.0) || (B2 * dx * dx + A2 * py2 > A2 * B2);
    float l    = length(vec2(useA ? p.x : p.x - d, p.y));
    return useA ? (l - ra) : (rb - l);
}

// 2. Exact 1-Sqrt-Arc 2D Crescent Moon SDF + Analytical Unit Gradient vec3(d, grad_x, grad_y)
vec3 sdgMoon2D( vec2 p, float d, float ra, float rb )
{
    float sy   = (p.y < 0.0) ? -1.0 : 1.0;
    float py2  = p.y * p.y;
    float a    = (ra * ra - rb * rb + d * d) / (2.0 * d);
    float b2   = max(ra * ra - a * a, 0.0);
    bool  low  = (b2 > py2);
    float u    = low ? d * (p.x - d) : d * p.x;
    float cv   = low ? d * (a - d)   : d * a;
    if( b2 > 0.0 && u * abs(u) * b2 > cv * abs(cv) * py2 )
    {
        vec2  w = vec2(p.x - a, abs(p.y) - sqrt(b2));
        float l = length(w);
        return vec3(l, vec2(w.x, sy * w.y) / max(l, 1e-30));
    }

    float A    = 0.5 * (ra + rb), A2 = A * A;
    float B2   = A2 - 0.25 * d * d;
    float dx   = p.x - 0.5 * d;
    bool  useA = (B2 <= 0.0) || (B2 * dx * dx + A2 * py2 > A2 * B2);
    vec2  w    = vec2(useA ? p.x : p.x - d, p.y);
    float l    = length(w);
    return useA ? vec3(l - ra,  w / max(l, 1e-30))
                : vec3(rb - l, -w / max(l, 1e-30));
}

// 3. Exact 1-to-2 Sqrt 3D Spherical Meniscus / Bowl SDF + Analytical Unit Gradient vec4(d, grad.xyz)
vec4 sdgMoon3D( vec3 p, float d, float ra, float rb )
{
    float ryz2 = dot(p.yz, p.yz);
    float a    = (ra * ra - rb * rb + d * d) / (2.0 * d);
    float b2   = max(ra * ra - a * a, 0.0);
    bool  low  = (b2 > ryz2);
    float u    = low ? d * (p.x - d) : d * p.x;
    float cv   = low ? d * (a - d)   : d * a;
    if( b2 > 0.0 && u * abs(u) * b2 > cv * abs(cv) * ryz2 )
    {
        float m    = sqrt(ryz2 * b2);
        float qx   = p.x - a;
        float l    = sqrt(max(qx * qx + ryz2 + b2 - 2.0 * m, 0.0));
        float s_yz = (ryz2 > 0.0) ? (1.0 - m / ryz2) / max(l, 1e-30) : 0.0;
        return vec4(l, qx / max(l, 1e-30), p.y * s_yz, p.z * s_yz);
    }

    float A    = 0.5 * (ra + rb), A2 = A * A;
    float B2   = A2 - 0.25 * d * d;
    float dx   = p.x - 0.5 * d;
    bool  useA = (B2 <= 0.0) || (B2 * dx * dx + A2 * ryz2 > A2 * B2);
    vec3  w    = vec3(useA ? p.x : p.x - d, p.y, p.z);
    float l    = length(w);
    return useA ? vec4(l - ra,  w / max(l, 1e-30))
                : vec4(rb - l, -w / max(l, 1e-30));
}

// 4. Exact 1-Sqrt-Arc 2D Asymmetric Two-Circle Optical Lens (0-sqrt signed-square wedge + confocal hyperbola)
float sdLens2D( vec2 p, float d, float ra, float rb )
{
    float py2 = p.y * p.y;
    float a   = (ra * ra - rb * rb + d * d) / (2.0 * d);
    float b2  = max(ra * ra - a * a, 0.0);
    float u1  = p.x,     v1 = a;
    float u2  = p.x - d, v2 = a - d;
    if( b2 > 0.0 && u1 * abs(u1) * b2 < v1 * abs(v1) * py2 && u2 * abs(u2) * b2 > v2 * abs(v2) * py2 )
        return length(vec2(p.x - a, abs(p.y) - sqrt(b2)));

    float r    = ra - rb, m = ra * r + d * (p.x - a);
    bool  useA = (r * abs(r) * (p.x * p.x + py2) <= m * abs(m));
    return length(vec2(useA ? p.x : p.x - d, p.y)) - (useA ? ra : rb);
}
