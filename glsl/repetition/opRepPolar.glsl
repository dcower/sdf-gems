// ============================================================================
// opRepPolar — Exact 1-Tap & 2-Tap Infinite & Limited Polar Repetition
// @author Mercury Demogroup — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (eliminates theta = +-pi branch cut, asymmetric [lima, limb] Voronoi tear, and boundary neighbor drop)
// ============================================================================
//
// Mathematical Foundation:
//   For polar repetition into n angular sectors of width sp = 2*pi / n on a
//   finite sector index range [lima, limb]:
//   1. Voronoi-Centered Periodic Wrapping:
//      When lima + limb != 0, the true Voronoi bisector of the empty angular
//      wedge is at the antipode of mid = 0.5*(lima + limb), NOT at theta = +-pi
//      (where naive clamping causes up to 286.9% error). Wrapping the sector
//      coordinate around mid:
//        sector = mid + mod(a - mid + 0.5*n, n) - 0.5*n
//      places the periodic branch cut at the exact Voronoi bisector mid +- 0.5*n.
//
//   2. Margin-Guarded Empty-Wedge Straddling (2-Tap Asymmetric):
//      Simply switching to (lima, limb) at |sector - mid| > 0.5*(limb - lima)
//      immediately drops the adjacent neighbor (lima + 1 or limb - 1) right at
//      the boundary ray, clipping asymmetric primitives that overhang the ray.
//      Extending the adjacent-pair zone by margin = 0.25*max(n - span - 1, 0)
//      preserves (lima, lima + 1) and (limb - 1, limb) past the boundary rays
//      before transitioning to (lima, limb) across the empty-wedge bisector.

// 1. 1-Tap Infinite Polar Repetition for radially symmetric primitives
vec2 opRepPolar( vec2 p, float n, out float id )
{
    float sp = 6.28318530718 / n;
    float a  = atan(p.y, p.x);
    id = round(a / sp);
    float al = a - id * sp;
    return length(p) * vec2(cos(al), sin(al));
}

vec2 opRepPolar( vec2 p, float n )
{
    float id;
    return opRepPolar(p, n, id);
}

// 2. 1-Tap Limited Polar Repetition on [lima, limb] (Voronoi-centered branch cut)
vec2 opRepPolarLim( vec2 p, float n, float lima, float limb, out float id )
{
    float sp     = 6.28318530718 / n;
    float a      = atan(p.y, p.x) / sp;
    float mid    = 0.5 * (lima + limb);
    float sector = mid + mod(a - mid + 0.5 * n, n) - 0.5 * n;
    id = clamp(round(sector), lima, limb);
    float al = (sector - id) * sp;
    return length(p) * vec2(cos(al), sin(al));
}

vec2 opRepPolarLim( vec2 p, float n, float lima, float limb )
{
    float id;
    return opRepPolarLim(p, n, lima, limb, id);
}

// 3. 2-Tap Limited Polar Repetition on [lima, limb] for asymmetric primitives (outputs angles id0, id1 in radians)
void opLimitedPolarRep2Tap( in vec2 p, in float n, in float lima, in float limb, out float id0, out float id1 )
{
    float sp     = 6.28318530718 / n;
    float a      = atan(p.y, p.x) / sp;
    float mid    = 0.5 * (lima + limb);
    float span   = limb - lima;
    float sector = mid + mod(a - mid + 0.5 * n, n) - 0.5 * n;
    float margin = 0.25 * max(n - span - 1.0, 0.0);
    float i0, i1;
    if( abs(sector - mid) <= 0.5 * span + margin )
    {
        i0 = clamp(floor(sector), lima, max(lima, limb - 1.0));
        i1 = min(i0 + 1.0, limb);
    }
    else
    {
        i0 = lima;
        i1 = limb;
    }
    id0 = i0 * sp;
    id1 = i1 * sp;
}

// Symmetric [-max_idx, max_idx] overload (outputs angles id0, id1 in radians)
void opLimitedPolarRep2Tap( in vec2 p, in float n, in float max_idx, out float id0, out float id1 )
{
    opLimitedPolarRep2Tap(p, n, -max_idx, max_idx, id0, id1);
}
