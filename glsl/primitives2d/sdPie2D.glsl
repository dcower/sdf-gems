// ============================================================================
// sdPie2D & sdgPie2D — Exact 0/1-Sqrt 2D Pie Sector & Analytical Unit Gradient
// @author Inigo Quilez — upgraded to 0/1-sqrt via Parabolic v^2-Cancellation using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 0 sqrts in exterior & interior radial-edge Voronoi region, 1 SFU elsewhere, 0 NaNs)
// ============================================================================

// 1. Exact 0/1-Sqrt 2D Circular Pie Sector:
//    By the Parabolic v^2-Cancellation Theorem, in sector coordinates u = dot(p, c),
//    v = c.y*|p.x| - c.x*p.y:
//    - In the exterior radial-edge normal strip (v >= 0 && 0 <= u <= r), distance is v (0 sqrts).
//    - In the interior (v < 0 && u >= 0), the straight radial edge (d = v) beats the circular arc
//      (d = sqrt(u^2 + v^2) - r) iff r + v > sqrt(u^2 + v^2) <=> r * (r + 2*v) > u^2 (0 sqrts!).
float sdPie2D( in vec2 p, in vec2 c, in float r )
{
    p.x = abs(p.x);
    float u = dot(p, c), v = c.y * p.x - c.x * p.y;
    if( u >= 0.0 && ((v >= 0.0) ? (u <= r) : (r + v > 0.0 && r * (r + 2.0 * v) > u * u)) )
        return v;
    float l = length((v > 0.0 && u > r) ? (p - r * c) : p);
    return (v > 0.0) ? l : ((u < 0.0 && l < 0.5 * r) ? -l : l - r);
}

// 2. Exact 0/1-SFU 2D Pie Sector Distance & Analytical Unit Gradient (down from 4 SFU ops in IQ's sdgPie, 0 NaNs)
vec3 sdgPie2D( in vec2 p, in vec2 c, in float r )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    float u = dot(p, c), v = c.y * p.x - c.x * p.y;
    if( u >= 0.0 && ((v >= 0.0) ? (u <= r) : (r + v > 0.0 && r * (r + 2.0 * v) > u * u)) )
        return vec3(v, sx * c.y, -c.x);
    vec2  q  = (v > 0.0 && u > r) ? (p - r * c) : p;
    float l2 = dot(q, q), il = inversesqrt(max(l2, 1e-30)), l = l2 * il;
    bool  na = (v <= 0.0 && u < 0.0 && l2 < 0.25 * r * r);
    float d  = (v > 0.0) ? l : (na ? -l : l - r);
    float sg = na ? -1.0 : 1.0;
    return vec3(d, sx * sg * q.x * il, sg * q.y * il);
}
