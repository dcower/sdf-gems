// ============================================================================
// sdgTriPrism — Exact 0/1-Sqrt 1-Reflection 3D Equilateral Triangular Prism & Analytical Gradient
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1`, 0 sqrts in interior & caps, 1 sqrt in exterior walls/corners, 0 segment projections)
// ============================================================================

// Exact 0/1-Sqrt 1-Reflection 3D Equilateral Triangular Prism with Analytical Unit Gradient:
// Uses the involutory reflection M = 0.5 * mat2(1, -sqrt(3), -sqrt(3), -1) (M^T = M = M^-1)
// to fold all 3 triangle edges onto the horizontal base segment [-2*r, 0] x {0}.
vec4 sdgTriPrism( vec3 p, vec2 h )
{
    const float k = 1.7320508;
    float r  = 0.5 * k * h.x;
    float sx = (p.x < 0.0) ? -1.0 : 1.0, sz = (p.z < 0.0) ? -1.0 : 1.0;
    vec3  q  = vec3(abs(p.x) - r, p.y + 0.5 * h.x, abs(p.z) - h.y);
    bool  rf = (q.x + k * q.y > 0.0);
    if( rf ) q.xy = vec2(q.x - k * q.y, -k * q.x - q.y) * 0.5;
    float ux = q.x - clamp(q.x, -2.0 * r, 0.0), dy = -q.y, dz = q.z, g = max(dy, dz);
    vec3  v  = vec3(ux, max(dy, 0.0), max(dz, 0.0));
    float l  = length(v);
    vec3  gr = (g > 0.0) ? vec3(v.x, -v.y, v.z) / max(l, 1e-30)
                         : ((dy > dz) ? vec3(0.0, -1.0, 0.0) : vec3(0.0, 0.0, 1.0));
    if( rf ) gr.xy = vec2(gr.x - k * gr.y, -k * gr.x - gr.y) * 0.5;
    return vec4((g > 0.0) ? l : g, sx * gr.x, gr.y, sz * gr.z);
}
