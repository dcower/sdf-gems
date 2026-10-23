// ============================================================================
// sdEgg2D, sdgEgg2D, sdEgg3D & sdgEgg3D — Exact 1-Length / 1-Sqrt-Cap Moss's Egg & 3D Egg of Revolution
// @author Inigo Quilez (sdEgg) — upgraded with confocal-hyperbola cap partition & signed-square cap tests using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 1 sqrt in caps & 2 in side arc down from 3 in 2D and 4 in 3D, `0.0` error)
// ============================================================================

// 1. Exact 2D Moss's Egg SDF (1 sqrt in bottom & top caps, 2 sqrts in side arc; 1 length() call total):
//    he = center separation, ra = bottom radius, rb = top radius, bu = bulge parameter in (0, 1].
float sdEgg2D( in vec2 p, in float he, in float ra, in float rb, in float bu )
{
    float r   = 0.5 * (he + ra + rb) / bu;
    float da  = r - ra;
    float db  = r - rb;
    float y   = (db * db - da * da - he * he) / (2.0 * he);
    float x2  = max(da * da - y * y, 0.0);

    p.x       = abs(p.x);
    float px2 = p.x * p.x;
    float pyb = p.y - he;
    float yb  = y + he;

    bool bot  = (p.y * abs(p.y) * x2 < y * abs(y) * px2);
    bool top  = (pyb * abs(pyb) * x2 > yb * abs(yb) * px2);

    vec2  q   = bot ? p : (top ? vec2(p.x, pyb) : p + vec2(sqrt(x2), y));
    float cr  = bot ? ra : (top ? rb : r);
    return length(q) - cr;
}

// 2. Exact 2D Moss's Egg SDF + Analytical Unit Gradient vec3(d, grad_x, grad_y)
vec3 sdgEgg2D( in vec2 p, in float he, in float ra, in float rb, in float bu )
{
    float r   = 0.5 * (he + ra + rb) / bu;
    float da  = r - ra;
    float db  = r - rb;
    float y   = (db * db - da * da - he * he) / (2.0 * he);
    float x2  = max(da * da - y * y, 0.0);

    float sx  = (p.x < 0.0) ? -1.0 : 1.0;
    p.x       = abs(p.x);
    float px2 = p.x * p.x;
    float pyb = p.y - he;
    float yb  = y + he;

    bool bot  = (p.y * abs(p.y) * x2 < y * abs(y) * px2);
    bool top  = (pyb * abs(pyb) * x2 > yb * abs(yb) * px2);

    vec2  q   = bot ? p : (top ? vec2(p.x, pyb) : p + vec2(sqrt(x2), y));
    float cr  = bot ? ra : (top ? rb : r);
    float l   = length(q);
    vec2  g   = (l > 1e-20) ? (vec2(sx * q.x, q.y) / l) : vec2(0.0, bot ? -1.0 : 1.0);
    return vec3(l - cr, g);
}

// 3. Exact 3D Moss's Egg Solid of Revolution + Analytical Unit Gradient vec4(d, grad.xyz):
//    Uses Radial-Toroidal Product Fusion m = sqrt(rxz2 * x2) (1 sqrt in caps, 2 sqrts in side arc, down from 4).
vec4 sdgEgg3D( in vec3 p, in float he, in float ra, in float rb, in float bu )
{
    float r    = 0.5 * (he + ra + rb) / bu;
    float da   = r - ra;
    float db   = r - rb;
    float y    = (db * db - da * da - he * he) / (2.0 * he);
    float x2   = max(da * da - y * y, 0.0);

    float rxz2 = dot(p.xz, p.xz);
    float pyb  = p.y - he;
    float yb   = y + he;

    if( p.y * abs(p.y) * x2 < y * abs(y) * rxz2 )
    {
        float l = sqrt(rxz2 + p.y * p.y);
        return vec4(l - ra, (l > 1e-20) ? (p / l) : vec3(0.0, -1.0, 0.0));
    }
    if( pyb * abs(pyb) * x2 > yb * abs(yb) * rxz2 )
    {
        float l = sqrt(rxz2 + pyb * pyb);
        return vec4(l - rb, (l > 1e-20) ? (vec3(p.x, pyb, p.z) / l) : vec3(0.0, 1.0, 0.0));
    }
    float m    = sqrt(rxz2 * x2);
    float qy   = p.y + y;
    float l    = sqrt(rxz2 + x2 + 2.0 * m + qy * qy);
    float s_xz = (rxz2 > 0.0) ? (1.0 + m / rxz2) / max(l, 1e-30) : 0.0;
    return vec4(l - r, p.x * s_xz, qy / max(l, 1e-30), p.z * s_xz);
}

// 4. Exact 3D Moss's Egg Solid of Revolution SDF
float sdEgg3D( in vec3 p, in float he, in float ra, in float rb, in float bu )
{
    return sdgEgg3D(p, he, ra, rb, bu).x;
}
