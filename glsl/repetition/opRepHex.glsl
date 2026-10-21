// ============================================================================
// opRepHex — 1-Tap Linear-Norm 2D Hexagonal (A2) Voronoi Domain Repetition
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (symmetric SDFs; 1 `round`, 0 `dot2` calls via |d.x| + 3|d.y| > 1.0)
// ============================================================================

vec2 opRepHex( vec2 p, float s )
{
    const vec2 k = vec2(1.0, 1.73205081);
    vec2 u = p / (s * k), i = round(u), d = u - i;
    if( abs(d.x) + 3.0 * abs(d.y) > 1.0 )
        i += vec2(d.x >= 0.0 ? 0.5 : -0.5, d.y >= 0.0 ? 0.5 : -0.5);
    return p - s * i * k;
}
