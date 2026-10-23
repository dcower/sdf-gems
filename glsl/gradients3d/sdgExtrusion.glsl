// ============================================================================
// sdgExtrusion — Exact 2D-to-3D Extrusion Operators with Analytical Gradients
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` across exterior bevels, flat caps, and interior)
// ============================================================================

#ifndef SDF_GEMS_SDG_EXTRUSION_GLSL
#define SDF_GEMS_SDG_EXTRUSION_GLSL

// 1. Sharp Linear Extrusion of any 2D SDF+gradient d2g = vec3(d, gx, gy)
vec4 sdgExtrusion( vec3 d2g, float z, float h )
{
    float u   = d2g.x, v = abs(z) - h;
    float sz  = (z < 0.0) ? -1.0 : 1.0;
    float g   = max(u, v);
    vec2  q   = max(vec2(u, v), 0.0);
    float l   = length(q);
    vec2  duv = (g > 0.0) ? q / l : ((u > v) ? vec2(1.0, 0.0) : vec2(0.0, 1.0));
    return vec4((g > 0.0) ? l : g, duv.x * d2g.yz, duv.y * sz);
}

// 2. Round-Beveled Extrusion with Analytical 3D Gradient
vec4 sdgExtrudeRound( vec3 d2g, float z, float h, float r )
{
    vec4 res = sdgExtrusion(vec3(d2g.x + r, d2g.yz), z, h - r);
    return vec4(res.x - r, res.yzw);
}

// 3. 45-Degree Chamfered Extrusion with Analytical 3D Gradient
vec4 sdgExtrudeChamfer( vec3 d2g, float z, float h, float c )
{
    float u   = d2g.x, v = abs(z) - h;
    float sz  = (z < 0.0) ? -1.0 : 1.0;
    float a   = min(u, v), b = max(u, v);
    const float s2 = 0.70710678;
    bool  w   = (b > 0.0 && b - a > c);
    float qa  = max(a + c, 0.0), l = length(vec2(qa, b));
    bool  ch  = (a + b + c) * s2 > b;
    float d   = w ? l : (ch ? (a + b + c) * s2 : b);
    vec2  ab  = w ? vec2(qa, b) / l : (ch ? vec2(s2) : vec2(0.0, 1.0));
    vec2  duv = (u < v) ? ab : ab.yx;
    return vec4(d, duv.x * d2g.yz, duv.y * sz);
}

// 4. n-Tier Ziggurat Stepped Extrusion with Analytical 3D Gradient
vec4 sdgExtrudeStairs( vec3 d2g, float z, float h, float r, float n )
{
    float u   = d2g.x, v = abs(z) - h;
    float sz  = (z < 0.0) ? -1.0 : 1.0;
    float s   = r / n, w = (v - u + r) / (2.0 * s);
    float j0  = clamp(floor(w), 0.0, n - 1.0) * s;
    float k1  = clamp(round(w), 1.0, n - 1.0) * s;
    vec2  p   = vec2(u + j0, v + r - s - j0);
    vec2  q   = vec2(-u - k1, -v - r + k1);
    bool  m   = max(p.x, p.y) < 0.0;
    vec2  c   = max(m ? q : p, 0.0);
    float l   = length(c), g = max(u, v);
    vec2  duv = (m && g > -l) ? ((u > v) ? vec2(1.0, 0.0) : vec2(0.0, 1.0))
                              : ((c.x > 0.0 && c.y > 0.0) ? (c / l) : (((m ? q.x : p.x) > (m ? q.y : p.y)) ? vec2(1.0, 0.0) : vec2(0.0, 1.0)));
    return vec4(m ? max(g, -l) : l, duv.x * d2g.yz, duv.y * sz);
}

// 5. Concave Cove / Cavetto Rabbet Extrusion with Analytical 3D Gradient
vec4 sdgExtrudeCove( vec3 d2g, float z, float h, float r )
{
    float u   = d2g.x, v = abs(z) - h;
    float sz  = (z < 0.0) ? -1.0 : 1.0;
    float a   = min(u, v), b = max(u, v);
    vec2  q   = (b > 0.0) ? vec2(max(a + r, 0.0), b) : vec2(a, b);
    float l   = length(q);
    float d   = (b > 0.0) ? l : max(b, r - l);
    vec2  ab  = (b > 0.0) ? q / l : ((r - l > b && l > 0.0) ? -q / l : vec2(0.0, 1.0));
    vec2  duv = (u < v) ? ab : ab.yx;
    return vec4(d, duv.x * d2g.yz, duv.y * sz);
}

// 6. Convex Rim Bead / Astragal Torus Extrusion with Analytical 3D Gradient
vec4 sdgExtrudeBead( vec3 d2g, float z, float h, float r )
{
    float u   = d2g.x, v = abs(z) - h;
    float sz  = (z < 0.0) ? -1.0 : 1.0;
    float a   = min(u, v), b = max(u, v);
    vec2  q   = (b < 0.0) ? vec2(max(a + r, 0.0), b) : vec2(a, b);
    float l   = length(q);
    bool  arc = (b >= 0.0) && (a > -r || l - r < b);
    float d   = (b < 0.0) ? -l : (arc ? l - r : b);
    vec2  ab  = (b < 0.0) ? -q / l : ((arc && l > 0.0) ? q / l : vec2(0.0, 1.0));
    vec2  duv = (u < v) ? ab : ab.yx;
    return vec4(d, duv.x * d2g.yz, duv.y * sz);
}

// 7. Tangent Roman Ogee (Cyma Recta) S-Curve Extrusion with Analytical 3D Gradient
vec4 sdgExtrudeOgee( vec3 d2g, float z, float h, float r )
{
    float u  = d2g.x, v = abs(z) - h, R = 0.5 * r;
    float sz = (z < 0.0) ? -1.0 : 1.0;
    float x  = u + R, y = v + R;
    if( x > 0.0 && y <= -R ) return vec4(u, d2g.yz, 0.0);
    bool  c  = (x <= 0.0) ? (y > R) : (y > 0.0 && 3.0 * y * y - 2.0 * y * (2.0 * x + R) > R * R);
    vec2  q  = vec2(c ? max(x + R, 0.0) : x, (c || x <= 0.0) ? y - R : y + R);
    float l  = length(q);
    if( c ) return vec4(l, (q.x / l) * d2g.yz, (q.y / l) * sz);
    if( x > 0.0 ) return vec4(l - R, (q.x / l) * d2g.yz, (q.y / l) * sz);
    bool  cv = (R - l > max(u, v)) && (l > 0.0);
    float d  = cv ? R - l : max(u, v);
    vec2  g  = cv ? -q / l : ((u > v) ? vec2(1.0, 0.0) : vec2(0.0, 1.0));
    return vec4(d, g.x * d2g.yz, g.y * sz);
}

#endif // SDF_GEMS_SDG_EXTRUSION_GLSL
