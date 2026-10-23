// ============================================================================
// opRevolutionLim — Exact 1-Tap Limited-Angle Partial Revolution of Any 2D SDF
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` across R^3 for all theta in (0, pi) in 1 tap of sdf2d)
// ============================================================================

// Step 1: Map 3D point p to 2D cross-section coordinates q.xy and cap coordinate q.z
// o  = major revolution radius (o > 0)
// sc = vec2(sin(theta), cos(theta)) for half-angle aperture theta in (0, pi)
vec3 opRevolutionLim( vec3 p, float o, vec2 sc )
{
    vec2  q = vec2(abs(p.x), p.z);
    float u = q.x * sc.x + q.y * sc.y;
    float w = q.x * sc.y - q.y * sc.x;
    float l = length(q);
    return vec3(((w > 0.0) ? u : l) - o, p.y, (w <= 0.0 && u < 0.0) ? -l : w);
}

// Step 2: Combine evaluated 2D SDF d2d = sdf2d(q.xy) with cap coordinate w = q.z (flat end-caps)
float opRevolutionLimFinalize( float d2d, float w )
{
    vec2 d = vec2(d2d, w);
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
}

// Optional Step 2b: Round-capped partial revolution for unsigned spine d_spine >= 0 (subtract radius r outside)
float opRevolutionArcFinalize( float d_spine, float w )
{
    return (w > 0.0) ? length(vec2(max(d_spine, 0.0), w)) + min(d_spine, 0.0) : d_spine;
}

// Overload for signed 2D SDF d2d with cap rounding radius r (||grad d|| = 1 everywhere)
float opRevolutionArcFinalize( float d2d, float w, float r )
{
    float ds = d2d + r;
    return ((w > 0.0) ? length(vec2(max(ds, 0.0), w)) + min(ds, 0.0) : ds) - r;
}
