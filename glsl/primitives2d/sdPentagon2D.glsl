// ============================================================================
// sdPentagon2D & sdPentagonalPrism — 0/1-Sqrt 2D Pentagon & 1-Sqrt 3D Prism
// @author Inigo Quilez (2D) — upgraded to 0/1-sqrt & 1-sqrt 3D prism using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in interior & 5 edge strips; `1` sqrt in 3D prism)
// ============================================================================

// 1. 0/1-Sqrt 2D Regular Pentagon (inradius r)
//    By the Fundamental-Chamber Vertex-Dominance Theorem, every point in the D5 chamber with
//    d.x = |p.x| - r*k.z > 0 has d.y = p.y - r > 0 (strictly exterior).
//    Eliminates clamp(), sign(), and the sqrt across 100% of the interior and all 5 edge strips!
float sdPentagon2D( in vec2 p, in float r )
{
    const vec3 k = vec3(0.809016994, 0.587785252, 0.726542528);
    p.x = abs(p.x);
    p -= 2.0 * min(dot(vec2(-k.x, k.y), p), 0.0) * vec2(-k.x, k.y);
    p -= 2.0 * min(dot(vec2( k.x, k.y), p), 0.0) * vec2( k.x, k.y);
    vec2 d = vec2(abs(p.x) - r * k.z, p.y - r);
    return (d.x > 0.0) ? length(d) : d.y;
}

// 2. Exact 1-Sqrt 3D Regular Pentagonal Prism (inradius h.x, half-height h.y)
float sdPentagonalPrism( in vec3 p, in vec2 h )
{
    const vec3 k = vec3(0.809016994, 0.587785252, 0.726542528);
    p.x = abs(p.x);
    p.xy -= 2.0 * min(dot(vec2(-k.x, k.y), p.xy), 0.0) * vec2(-k.x, k.y);
    p.xy -= 2.0 * min(dot(vec2( k.x, k.y), p.xy), 0.0) * vec2( k.x, k.y);
    vec3 d = vec3(abs(p.x) - h.x * k.z, p.y - h.x, abs(p.z) - h.y);
    return min(max(d.y, d.z), 0.0) + length(max(d, 0.0));
}
