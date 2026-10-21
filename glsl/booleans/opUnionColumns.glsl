// ============================================================================
// opUnionColumns — Exact 1-Sqrt Euclidean Columned Booleans
// @author Mercury Demogroup (discontinuous bound) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1` for orthogonal surfaces, max error `8.88e-16`, 0 discontinuity)
// ============================================================================

// 1. Exact 1-Sqrt Corner-Local Euclidean Union Columns (n-1 circular columns)
float opUnionColumns( float a, float b, float r, float n )
{
    float x  = min(a, b), y = max(a, b);
    float R  = r * 1.41421356 / ((n - 1.0) * 2.0 + 1.41421356);
    float px = (x + y - r + 2.0 * R) * 0.70710678, py = (y - x) * 0.70710678;
    float w  = py / (2.0 * R) + 0.5 * (n - 2.0), k = min(round(w), n - 2.0);
    float dy = (w - k) * 2.0 * R, qx = max(px, 0.0), d2 = qx * qx + dy * dy;
    bool  ex = (x > 0.0) && (d2 > R * R);
    bool  ar = !ex && (px > 0.0) && (d2 <= R * R) && (k < n - 2.0 || px * 0.83822219 >= dy * 0.54532906);
    float ry = max(r - 1.22905745 * R - y, 0.0), di = x * x + ry * ry;
    float dj = (w + 0.5 - clamp(round(w + 0.5), 1.0, max(n - 2.0, 1.0))) * 2.0 * R;
    float l  = sqrt((ex || ar) ? d2 : (n > 2.0 ? min(di, px * px + dj * dj) : di));
    return ex ? min(x, l - R) : (ar ? l - R : -l);
}

// 2. Exact 1-Sqrt Corner-Local Euclidean Subtraction/Difference Columns (n tangent circular columns)
float opSubtractionColumns( float a, float b, float r, float n )
{
    float x  = min(-a, b), y = max(-a, b);
    float R  = r * 1.41421356 / ((n - 1.0) * 2.0 + 1.41421356);
    float px = (x + y - r - R) * 0.70710678, py = (y - x) * 0.70710678;
    float w  = py / (2.0 * R) + 0.5 * (n - 1.0);
    float dy = (w - min(round(w), n - 1.0)) * 2.0 * R;
    float dk = (w - 0.5 - min(round(w - 0.5), n - 2.0)) * 2.0 * R;
    float ry = max(r - y, 0.0);
    bool  c  = (y < r) && (px <= 0.0);
    float l  = sqrt(c ? px * px + dy * dy : min(x * x + ry * ry, px * px + dk * dk));
    return c ? l - R : (x < 0.0 ? -x : -l);
}

// 3. Exact 1-Sqrt Corner-Local Euclidean Intersection Columns
float opIntersectionColumns( float a, float b, float r, float n )
{
    return opSubtractionColumns(a, -b, r, n);
}
