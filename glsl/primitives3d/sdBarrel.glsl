// ============================================================================
// sdBarrel — Exact 3D Convex Barrel & Concave Waisted Spool
// @author Dillon Cower (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdBarrel( vec3 p, float r1, float r2, float h )
{
    float b  = r2 - r1;
    float sb = sign(b);
    float d  = 0.5*(h*h - b*b)/b;
    float R  = d + b;
    vec2  q  = vec2( length(p.xz), abs(p.y) );
    float u  = q.x - (r1 - d);
    float s  = max( sb*(d*q.y - h*u), sb*((R - h)*(R - h + 2.0*q.y) - u*u) );
    return (s < 0.0) ? sb*length(vec2(u, q.y)) - R :
           (q.x < r1) ? q.y - h :
                        length(q - vec2(r1, h));
}
