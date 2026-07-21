/*  zzz live loop — COMMON tab (shared by Buffer A and Image).
 *  Self-paving primes<->zeros bootstrap, transliterated from loop.c /
 *  loop_frame.c (verified flawless to 615 primes in genuine fp32).
 *
 *  STATE LAYOUT in Buffer A (one RGBA32F texture, indexed linearly):
 *    linear index i  <->  texel ivec2(i % STRIDE, i / STRIDE)
 *    .r channel, index k          : Z[k] = gamma_{k+1}   (located zero)
 *    .g channel, index x (integer): 1.0 if x is prime    (indicator field)
 *    texel 0 .b/.a                : nz , Xknown           (control 0)
 *    texel 1 .b/.a                : kmin , iter           (control 1)
 *  Zeros and the prime-indicator share the grid but live in different
 *  channels, so there is no packing/compaction step.
 */
const float PI   = 3.14159265358979;
const float TWO_PI = 6.28318530717959;
const float EG   = 0.5772156649015329;
const float DTAU = 0.30;     // detector strong-peak threshold
const float DLO  = 0.12;     // ambiguous-band floor
const int   STRIDE = 256;    // texels per row of state
const int   ZBASE  = 0;      // zeros start at linear index 0 (.r channel)
const int   NZMAX  = 60000;  // zero-list cap
const int   NSEED  = 40;
const float KFIX   = 5.0;    // kappa margin (fixed; over-supplied, robust)

// first 40 nontrivial ordinates (the finite seed; not part of the claim)
float seed(int i){
    float s[40] = float[](
    14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073,
    48.005151,49.773832,52.970321,56.446248,59.347044,60.831779,65.112544,67.079811,
    69.546402,72.067158,75.704691,77.144840,79.337375,82.910381,84.735493,87.425275,
    88.809111,92.491899,94.651344,95.870634,98.831194,101.317851,103.725538,105.446623,
    107.168611,111.029536,111.874659,114.320221,116.226680,118.790783,121.370125,122.946829);
    return s[i];
}

ivec2 idx2(int i){ return ivec2(i % STRIDE, i / STRIDE); }

// fp32 phase-safe cos / sin (Cody-Waite 2pi split)
float pcos(float a){ float k=floor(a*0.15915494+0.5); float r=(a-k*6.2831855)+k*1.7484555e-7; return cos(r); }
float psin(float a){ return pcos(a-1.5707963); }

float N0f(float t){ float u=t/TWO_PI; return u*log(u/2.718281828)+0.875; }
float Eif(float z){ float s=EG+log(z), t=1.0; for(int k=1;k<64;k++){ t*=z/float(k); s+=t/float(k);} return s; }
float envf(float x, int n){ return (20.0/4.0)*Eif(log(float(n)-4.0))*log(x)/x; } // (20/4)Li(n-4) log x / x

// invert N0(t)=k-1/2  (smooth gamma estimate; seeds the parallel bisection)
float smoothGamma(int k){
    float tgt=float(k)-0.5, t=TWO_PI*tgt/log(tgt+2.0)+10.0;
    for(int i=0;i<40;i++){ float f=N0f(t)-tgt, d=log(t/TWO_PI)/TWO_PI; if(d<1e-9)d=1e-9; t-=f/d; if(t<1.0)t=1.0; }
    return t;
}

// max zero index resolvable from primes<=X at margin kappa (n_frontier)
int nFrontier(int X, float kmin){
    if(X<2) return 0;
    float T=TWO_PI*pow(float(X), TWO_PI/kmin), u=T/TWO_PI;
    float n=floor(u*(log(u)-1.0));
    return (n>float(NZMAX-2))? NZMAX-2 : int(n);
}

const float HOLD = 1.2;     // seconds each prime is held on screen

// current prime + scroll window. The Common tab can't name pass uniforms, so the
// caller passes them in:  primeWindow(iChannel0, iTime, pc, xLo, xHi).
void primeWindow(sampler2D buf, float t, out int pc, out float xLo, out float xHi){
    vec4 c0 = texelFetch(buf, idx2(0),0);
    int Xk = int(c0.a+0.5); if(Xk<3) Xk=3;
    float nf = t/HOLD + 1.0;                     // CONTINUOUS (no floor => no skips)
    float Lg = log(nf+1.0);
    float pest = nf*(Lg + log(Lg+1.0));          // smooth, monotone, ~ p_n
    int cur = int(pest+0.5); if(cur>Xk) cur=Xk; if(cur<2) cur=2;
    pc=2; for(int d=0; d<80; d++){ int q=cur-d;
        if(q>=2 && texelFetch(buf, idx2(q),0).g>0.5){ pc=q; break; } }
    float pcf=float(pc);
    float WIN = max(36.0, 9.0*log(pcf+2.0));
    xHi = pcf + 0.22*WIN; xLo = xHi - WIN;
    if(xLo<2.0){ xLo=2.0; xHi=xLo+WIN; }
}
