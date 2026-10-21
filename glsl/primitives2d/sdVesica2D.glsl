// ============================================================================
// sdVesica2D — Exact 1-to-2 Sqrt 2D Vesica Piscis (Signed-Square Cone Reduction)
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt` in main body & 2 in tip cone, `0.0` error)
// ============================================================================

float sdVesica2D( vec2 p, float r, float d )
{
    p = abs(p);
    float b2 = r * r - d * d;
    float x  = p.x + d;
    return (d * d * p.y * p.y > b2 * x * x)
         ? length(vec2(p.x, p.y - sqrt(b2)))
         : length(vec2(x, p.y)) - r;
}

// Physical half-length (r) & half-width (w) parameterization: 1 `sqrt` everywhere in R^2!
float sdVesicaWidth2D( vec2 p, float r, float w )
{
    p = abs(p);
    float d = 0.5 * (r * r - w * w) / w;
    return (r * p.x < d * (p.y - r))
         ? length(vec2(p.x, p.y - r))
         : length(vec2(p.x + d, p.y)) - (d + w);
}
