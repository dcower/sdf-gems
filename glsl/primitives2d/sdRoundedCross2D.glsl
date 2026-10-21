// ============================================================================
// sdRoundedCross2D — Division-Free 2D Rounded Cross & 1-Sqrt 3D Rounded Cross Prism
// @author Inigo Quilez — upgraded via Division-Free Region Test & 1-Sqrt 3D Lifting using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1` everywhere; 0 divisions in region test, 1 dot2 in corner branch, 0/1-sqrt in 3D prism)
// ============================================================================
//
// Mathematical Foundation:
//   In Inigo Quilez's baseline `sdRoundedCross(vec2 p, float h)`
//   (`distfunctions2d.txt` lines 1334-1372):
//     float k = 0.5*(h + 1.0/h);
//     p = abs(p);
//     return (p.x < 1.0 && p.y < p.x*(k - h) + h)
//          ? k - sqrt(dot2(p - vec2(1, k)))
//          : sqrt(min(dot2(p - vec2(0, h)), dot2(p - vec2(1, 0))));
//
//   1. Division-Free Region Test & Single-Dot2 Corner Selector:
//      Since k - h = (1 - h^2) / (2*h) and h > 0, the region test
//        p.y < p.x * (k - h) + h
//      is algebraically equivalent to the division-free polynomial inequality:
//        2.0 * h * (p.y - h) < p.x * (1.0 - h * h),
//      deferring the reciprocal `1.0 / h` strictly to the concave arc branch!
//      Furthermore, in the exterior corner branch (`else`), comparing squared
//      distances to (0, h) and (1, 0):
//        ||p - (0, h)||^2 < ||p - (1, 0)||^2   <==>   2*p.x - 2*h*p.y < 1 - h^2,
//      selects the active corner with a single linear test (saving 1 dot2).
//
//   2. 1-Sqrt 3D Extruded `sdRoundedCrossPrism(vec3 p, float h, float hz)`:
//      Because the `else` branch of `sdRoundedCross2D` lies 100% in the 2D
//      exterior (d_2D >= 0) and its squared distance is a pure polynomial
//      dot2(q - c), the 3D extruded prism evaluates:
//        sqrt(dot2(q - c) + max(dz, 0.0)^2)
//      with 1 sqrt (instead of 2!) across the entire exterior corner region,
//      0 sqrts above/below the 2D interior (`dz > 0 && r2 >= k*k`), and 1 sqrt
//      inside the z-slab (`dz <= 0`).

// 1. Division-Free Region-Tested 2D Rounded Cross (1 sqrt, 1 dot2 in corner branch)
float sdRoundedCross2D( in vec2 p, in float h )
{
    p = abs(p);
    if( p.x < 1.0 && 2.0 * h * (p.y - h) < p.x * (1.0 - h * h) )
    {
        float k = 0.5 * (h + 1.0 / h);
        return k - length(p - vec2(1.0, k));
    }
    vec2 c = (2.0 * p.x - 2.0 * h * p.y < 1.0 - h * h) ? vec2(0.0, h) : vec2(1.0, 0.0);
    return length(p - c);
}

// 2. Analytical 2D Distance and Unit Gradient vec3(d, grad.xy) (1 inversesqrt everywhere)
vec3 sdgRoundedCross2D( in vec2 p, in float h )
{
    vec2 s = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    if( p.x < 1.0 && 2.0 * h * (p.y - h) < p.x * (1.0 - h * h) )
    {
        float k    = 0.5 * (h + 1.0 / h);
        vec2  v    = vec2(1.0, k) - p;
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        return vec3(k - v2 * invV, s * v * invV);
    }
    vec2  c    = (2.0 * p.x - 2.0 * h * p.y < 1.0 - h * h) ? vec2(0.0, h) : vec2(1.0, 0.0);
    vec2  v    = p - c;
    float v2   = dot(v, v);
    float invV = inversesqrt(max(v2, 1e-30));
    return vec3(v2 * invV, s * v * invV);
}

// 3. Exact 0/1-Sqrt 3D Rounded Cross Prism (0 sqrts above/below interior, 1 sqrt in z-slab & exterior corners)
float sdRoundedCrossPrism( in vec3 p, in float h, in float hz )
{
    vec2  q  = abs(p.xy);
    float dz = abs(p.z) - hz;
    if( q.x < 1.0 && 2.0 * h * (q.y - h) < q.x * (1.0 - h * h) )
    {
        float k  = 0.5 * (h + 1.0 / h);
        float r2 = dot(q - vec2(1.0, k), q - vec2(1.0, k));
        if( dz > 0.0 && r2 >= k * k ) return dz; // 0 sqrts above/below 2D interior!
        float d2 = k - sqrt(r2);
        return (dz <= 0.0) ? max(d2, dz) : length(vec2(d2, dz));
    }
    vec2  c  = (2.0 * q.x - 2.0 * h * q.y < 1.0 - h * h) ? vec2(0.0, h) : vec2(1.0, 0.0);
    float zc = max(dz, 0.0);
    return sqrt(dot(q - c, q - c) + zc * zc); // 1 sqrt!
}
