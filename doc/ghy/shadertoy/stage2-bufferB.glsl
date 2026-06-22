/*  zzz live loop — BUFFER B tab (trace cache).  Set its iChannel0 = Buffer A.
 *
 *  Computes the detector r(x) ONCE per screen column (row 0 only), so the Image
 *  pass can read the trace with a texelFetch instead of a 384-term cos sum per
 *  pixel. Uses the shared primeWindow() (Common) so its column->x mapping is
 *  identical to the Image pass.
 */
const int NVIS = 384;     // zeros used for the background trace (visual only)

float RdetVis(float x, int nz){
    int n = (nz<NVIS)? nz : NVIS;
    float lx=log(x), s=0.0;
    for(int k=0;k<NVIS;k++){ if(k>=n) break; s += pcos(texelFetch(iChannel0, idx2(k),0).r * lx); }
    float psi = 1.0 + 1.0/(x - x*x*x) - (4.0/sqrt(x))*s;
    return psi/envf(x, n);
}

void mainImage(out vec4 O, in vec2 F){
    if(int(F.y) != 0){ O = vec4(0.0); return; }          // only row 0 holds the trace
    int pc; float xLo, xHi; primeWindow(iChannel0, iTime, pc, xLo, xHi);
    int nz = int(texelFetch(iChannel0, idx2(0),0).b + 0.5);
    float x = xLo + (xHi-xLo) * (F.x / iResolution.x);   // this column's number-line position
    O = vec4(RdetVis(x, nz), 0.0, 0.0, 1.0);             // r(x) cached in .r
}
