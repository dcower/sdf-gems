// ============================================================================
// sdEquilateralTriangle2D — Exact 0/1-Sqrt 2D Equilateral Triangle
// @author Inigo Quilez — upgraded to 0-sqrt interior/face-strip using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in interior & face strips, `1` in corner wedges)
// ============================================================================
//
// Mathematical Derivation (Fundamental-Chamber Double Inequality):
//   After folding p via `p.x = abs(p.x) - r; p.y += r/k; if (p.x + k*p.y > 0) ...`,
//   the fundamental chamber is bounded by the two bisector rays p.x = -2*r + k*p.y
//   and p.x = -k*p.y, so every point in the chamber satisfies:
//     -2*r + k*p.y <= p.x <= -k*p.y   <==>   |p.x + r| - r <= -k * p.y.
//   Therefore, defining dx = abs(p.x + r) - r:
//     dx > 0  ==>  -k * p.y >= dx > 0  ==>  p.y < 0 (strictly in the 3 exterior corner wedges!).
//   Whenever dx <= 0 (the entire interior and all 3 exterior flat-edge normal strips),
//   the closest boundary point is (p.x, 0) at exact signed distance -p.y, eliminating
//   `clamp`, `sign`, and the `sqrt` (`length`) across the interior and face strips!

#ifndef SDF_GEMS_SD_EQUILATERAL_TRIANGLE2D_GLSL
#define SDF_GEMS_SD_EQUILATERAL_TRIANGLE2D_GLSL

float sdEquilateralTriangle2D( vec2 p, float r )
{
    const float k = sqrt(3.0);
    p.x = abs(p.x) - r;
    p.y = p.y + r/k;
    if( p.x + k*p.y > 0.0 ) p = vec2(p.x - k*p.y, -k*p.x - p.y)/2.0;
    float dx = abs(p.x + r) - r;
    return (dx > 0.0) ? length(vec2(dx, p.y)) : -p.y;
}

#endif // SDF_GEMS_SD_EQUILATERAL_TRIANGLE2D_GLSL
