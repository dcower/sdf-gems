// ============================================================================
// opGroove — Exact 1-Sqrt Euclidean Carpenter's Tongue & Groove Joint Booleans
// @author Mercury Demogroup (L-infinity bound) — upgraded to exact Euclidean & oblique Gram-Schmidt using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1`, `0.0` error vs `29.3%` orthogonal & `0.500` oblique error in hg_sdf)
// ============================================================================

// 1. Exact 4-Line 1-Sqrt Euclidean Carpenter's Tongue
float opTongue( float a, float b, float ra, float rb )
{
    vec2  q = vec2(a - ra, abs(b) - rb), i = vec2(-a, -q.y);
    bool  e = min(a, max(q.x, q.y)) > 0.0;
    float l = length(max(e ? q : i, 0.0));
    return e ? min(a, l) : -min(ra - a, l);
}

// 2. Exact 1-Sqrt Euclidean Carpenter's Groove (exact De Morgan dual of opTongue)
float opGroove( float a, float b, float ra, float rb )
{
    return -opTongue(-a, b, ra, rb);
}

// 3. Angle-Exact Oblique Carpenter's Tongue & Groove via Gram-Schmidt Surface-Anchored Orthonormal Map
// c = dot(grad(a), grad(b)) = cos(theta)
float opTongueOblique( float a, float b, float ra, float rb, float c )
{
    return opTongue(a, (b - c * a) * inversesqrt(max(1.0 - c * c, 1e-5)), ra, rb);
}

float opGrooveOblique( float a, float b, float ra, float rb, float c )
{
    return -opTongueOblique(-a, b, ra, rb, -c);
}
