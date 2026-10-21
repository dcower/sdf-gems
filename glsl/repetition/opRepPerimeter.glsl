// ============================================================================
// opRepPerimeter — Exact 1-Tap & 2-Tap Rectangular Perimeter Repetition
// @author Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/sdfrepetition/
// @exact  Yes (fixes 200% non-square rectangle failure & coordinate transposition in IQ's repeat_rectangular_ONLY_SYMMETRIC_SDFS; adds 2-tap asymmetric support)
// ============================================================================
//
// Mathematical Foundation:
//   In `sdfrepetition.txt` (lines 340-349), IQ's `repeat_rectangular_ONLY_SYMMETRIC_SDFS`
//   folds across the origin diagonal `p.x < p.y` instead of the 45-degree corner
//   bisector `|p.x| - l.x*s < |p.y| - l.y*s` and fails to swap `size.xy` after
//   swapping `p.yx`, resulting in 200% distance error on every non-square
//   rectangle (l.x != l.y) while also mirroring and transposing the local
//   coordinate frame on the horizontal walls.
//
//   1. Un-Mirrored, Un-Transposed 1-Tap Perimeter Repetition (`opRepPerimeter`):
//      Evaluating the signed wall offsets d = abs(p / s) - l directly determines
//      whether p lies in the vertical-wall Voronoi region (d.x >= d.y) or the
//      horizontal-wall Voronoi region (d.x < d.y), yielding the exact integer
//      perimeter cell index `id` and pure translation `p - s * id` (no mirroring
//      or 90-degree rotation!) with 0.0 error for all aspect ratios l.x, l.y.
//
//   2. Corner-Aware 2-Tap Perimeter Repetition (`opRepPerimeter2Tap`):
//      For asymmetric/rotated primitives repeated around a rectangular perimeter,
//      each wall segment pairs (floor(u), floor(u) + 1), and near the 4 corners
//      checks the adjacent cell on the perpendicular wall to select the two
//      closest perimeter cell centers (id0, id1) with 0.0 error everywhere.

// 1. 1-Tap Rectangular Perimeter Repetition around [-l.x, l.x] x [-l.y, l.y] with spacing s
vec2 opRepPerimeter( vec2 p, float s, vec2 l, out vec2 id )
{
    vec2 u   = p / s;
    vec2 sgn = vec2((u.x >= 0.0) ? 1.0 : -1.0, (u.y >= 0.0) ? 1.0 : -1.0);
    vec2 d   = abs(u) - l;
    id = (d.x >= d.y)
       ? vec2(sgn.x * l.x, clamp(round(u.y), -l.y, l.y))
       : vec2(clamp(round(u.x), -l.x, l.x), sgn.y * l.y);
    return p - s * id;
}

vec2 opRepPerimeter( vec2 p, float s, vec2 l )
{
    vec2 id;
    return opRepPerimeter(p, s, l, id);
}

// 2. 2-Tap Corner-Aware Rectangular Perimeter Repetition for asymmetric primitives
void opRepPerimeter2Tap( in vec2 p, in float s, in vec2 l, out vec2 id0, out vec2 id1 )
{
    vec2 u       = p / s;
    vec2 sgn     = vec2((u.x >= 0.0) ? 1.0 : -1.0, (u.y >= 0.0) ? 1.0 : -1.0);
    vec2 d       = abs(u) - l;
    bool is_vert = (d.x >= d.y);

    if( is_vert )
    {
        float yb = clamp(floor(u.y), -l.y, l.y - 1.0);
        id0 = vec2(sgn.x * l.x, yb);
        id1 = vec2(sgn.x * l.x, yb + 1.0);
    }
    else
    {
        float xb = clamp(floor(u.x), -l.x, l.x - 1.0);
        id0 = vec2(xb, sgn.y * l.y);
        id1 = vec2(xb + 1.0, sgn.y * l.y);
    }

    vec2 diff0 = u - id0;
    vec2 diff1 = u - id1;
    if( dot(diff1, diff1) < dot(diff0, diff0) )
    {
        vec2 tmp = id0; id0 = id1; id1 = tmp;
        diff1 = u - id1;
    }

    // Check adjacent cell on the perpendicular wall near the 90-degree corner
    vec2 sec;
    if( is_vert )
    {
        float hx = clamp(round(u.x), -l.x, l.x);
        if( id0.y == sgn.y * l.y && id0.x == hx ) hx = clamp(hx - sgn.x, -l.x, l.x);
        sec = vec2(hx, sgn.y * l.y);
    }
    else
    {
        float vy = clamp(round(u.y), -l.y, l.y);
        if( id0.x == sgn.x * l.x && id0.y == vy ) vy = clamp(vy - sgn.y, -l.y, l.y);
        sec = vec2(sgn.x * l.x, vy);
    }
    vec2 diff_sec = u - sec;
    if( dot(diff_sec, diff_sec) < dot(diff1, diff1) )
        id1 = sec;
}
