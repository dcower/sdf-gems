// ============================================================================
// opExtrudeMolding — Exact 1-Sqrt Cove, Bead & Roman Ogee Molding Extrusions
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt` total, `0.0` to `4.44e-16` max error)
// ============================================================================

#ifndef SDF_GEMS_OP_EXTRUDE_MOLDING_GLSL
#define SDF_GEMS_OP_EXTRUDE_MOLDING_GLSL

// 1. Exact 4-Line 1-Sqrt Quarter-Circle Concave Cove / Cavetto Molding Extrusion
float opExtrudeCove( in float d2d, in float z, in float h, in float r )
{
    float u = d2d, v = abs(z) - h;
    float a = min(u, v), b = max(u, v);
    float l = length((b > 0.0) ? vec2(max(a + r, 0.0), b) : vec2(a, b));
    return (b > 0.0) ? l : max(b, r - l);
}

// 2. Exact 4-Line 1-Sqrt Convex Rim Bead / Astragal Torus Molding Extrusion
//    Fixes the 29.3% interior quirk-corner error of naive min(d_extrude, length(vec2(u,v)) - r).
float opExtrudeBead( in float d2d, in float z, in float h, in float r )
{
    float u = d2d, v = abs(z) - h;
    float a = min(u, v), b = max(u, v);
    float l = length((b < 0.0) ? vec2(max(a + r, 0.0), b) : vec2(a, b));
    return (b < 0.0) ? -l : ((a > -r) ? l - r : min(b, l - r));
}

// 3. Exact 0/1-Sqrt Tangent Roman Ogee (Cyma Recta) S-Curve Molding Extrusion
//    Uses the Quadratic Apollonius Separator (3*y^2 - 2*y*(2*x + R) > R^2) to select
//    between the top-cap vertex and lower convex arc BEFORE taking the square root.
float opExtrudeOgee( in float d2d, in float z, in float h, in float r )
{
    float u = d2d, v = abs(z) - h, R = 0.5 * r;
    float x = u + R, y = v + R;
    if( x > 0.0 && y <= -R ) return u;
    bool  c = (x <= 0.0) ? (y > R) : (y > 0.0 && 3.0 * y * y - 2.0 * y * (2.0 * x + R) > R * R);
    vec2  q = vec2(c ? max(x + R, 0.0) : x, (c || x <= 0.0) ? y - R : y + R);
    float l = length(q);
    return c ? l : ((x <= 0.0) ? max(max(u, v), R - l) : l - R);
}

#endif // SDF_GEMS_OP_EXTRUDE_MOLDING_GLSL
