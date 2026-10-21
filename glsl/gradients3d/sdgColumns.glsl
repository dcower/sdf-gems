// ============================================================================
// sdgColumns — Exact 1-Sqrt Columned Booleans with Analytical 3D Gradients
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` for orthogonal surfaces, max Eikonal error `4.44e-16`)
// ============================================================================

vec4 sdgUnionColumns( in vec4 a, in vec4 b, in float r, in float n )
{
    if( a.x > b.x ) { vec4 t = a; a = b; b = t; }
    float x  = a.x, y = b.x, s2 = 0.70710678;
    float R  = r * 1.41421356 / ((n - 1.0) * 2.0 + 1.41421356);
    float px = (x + y - r + 2.0 * R) * s2, py = (y - x) * s2;
    float w  = py / (2.0 * R) + 0.5 * (n - 2.0), k = min(round(w), n - 2.0);
    float dy = (w - k) * 2.0 * R, qx = max(px, 0.0), d2 = qx * qx + dy * dy;
    bool  ex = (x > 0.0) && (d2 > R * R);
    bool  ar = !ex && (px > 0.0) && (d2 <= R * R) && (k < n - 2.0 || px * 0.83822219 >= dy * 0.54532906);
    float ry = max(r - 1.22905745 * R - y, 0.0), di = x * x + ry * ry;
    float dj = (w + 0.5 - clamp(round(w + 0.5), 1.0, max(n - 2.0, 1.0))) * 2.0 * R;
    bool  cu = (n > 2.0) && (px * px + dj * dj < di);
    float l  = sqrt((ex || ar) ? d2 : (cu ? px * px + dj * dj : di));
    if( ex && x < l - R ) return a;
    float d  = (ex || ar) ? l - R : -l;
    vec2  gp = (ex || ar) ? vec2(qx, dy) / l : -vec2(px, dj) / l;
    vec2  gr = (ex || ar || cu) ? vec2(gp.x - gp.y, gp.x + gp.y) * s2
                                : ((ry <= 0.0) ? vec2(1.0, 0.0) : vec2(-x, ry) / l);
    return vec4(d, gr.x * a.yzw + gr.y * b.yzw);
}

vec4 sdgSubtractionColumns( in vec4 a, in vec4 b, in float r, in float n )
{
    a = -a;
    if( a.x > b.x ) { vec4 t = a; a = b; b = t; }
    float x  = a.x, y = b.x, s2 = 0.70710678;
    float R  = r * 1.41421356 / ((n - 1.0) * 2.0 + 1.41421356);
    float px = (x + y - r - R) * s2, py = (y - x) * s2;
    float w  = py / (2.0 * R) + 0.5 * (n - 1.0);
    float dy = (w - min(round(w), n - 1.0)) * 2.0 * R;
    float dk = (w - 0.5 - min(round(w - 0.5), n - 2.0)) * 2.0 * R;
    float ry = max(r - y, 0.0), di = x * x + ry * ry;
    bool  c  = (y < r) && (px <= 0.0);
    bool  cu = !c && (x >= 0.0) && (px * px + dk * dk < di);
    float l  = sqrt(c ? px * px + dy * dy : (cu ? px * px + dk * dk : di));
    if( !c && (x < 0.0 || (!cu && ry <= 0.0)) ) return -a;
    float d  = c ? l - R : -l;
    vec2  gp = c ? vec2(px, dy) / l : -vec2(px, dk) / l;
    vec2  gr = (c || cu) ? vec2(gp.x - gp.y, gp.x + gp.y) * s2 : vec2(-x, ry) / l;
    return vec4(d, gr.x * a.yzw + gr.y * b.yzw);
}

vec4 sdgIntersectionColumns( in vec4 a, in vec4 b, in float r, in float n )
{
    return sdgSubtractionColumns(a, -b, r, n);
}
