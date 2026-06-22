/*  zzz live loop — BUFFER A tab.  Set this buffer's iChannel0 = Buffer A
 *  (feedback).  One mainImage call per texel per frame advances the loop
 *  exactly like one iteration of loop_frame.c: backward (extend the prime
 *  frontier from current zeros) + forward (locate a batch of new zeros by
 *  method-B bisection over the discovered primes).  See the Common tab for
 *  the state layout.  DEMO CAPS keep per-fragment work inside browser limits.
 */
const int BATCH = 12;     // zeros located per frame
const int BAND  = 4;      // prime frontier advanced per frame
const int NLIM  = 8000;   // demo zero cap (=> primes to ~3000); raise carefully
const int XKMAX = 4096;   // FB prime-scan bound

bool isPP(int x){                                   // perfect power p^m, m>=2 ?
    for(int b=2;b<14;b++){
        float a=floor(pow(float(x),1.0/float(b))+0.5);
        if(a<2.0) break;
        float p=1.0; for(int e=0;e<b;e++) p*=a;
        if(int(p+0.5)==x) return true;
    }
    return false;
}

// detector r(x) from the nz committed zeros (.r channel of Buffer A)
float detFire(float x, int nz){
    float lx=log(x), s=0.0;
    for(int k=0;k<NLIM;k++){ if(k>=nz) break; s += pcos(texelFetch(iChannel0, idx2(k),0).r * lx); }
    float psi = 1.0 + 1.0/(x - x*x*x) - (4.0/sqrt(x))*s;
    return psi/envf(x, nz);
}

// F_B(t): method-B count from the prime-indicator field (.g), primes <= Xk
float FB(float t, int Xk){
    float acc=0.0;
    for(int x=2;x<XKMAX;x++){
        if(x>Xk) break;
        if(texelFetch(iChannel0, idx2(x),0).g > 0.5){       // x is prime
            float lp=log(float(x)), pm=float(x);
            for(int m=1;m<20;m++){
                float hf=pow(float(x), -0.5*float(m));
                acc += (hf/float(m))*psin(float(m)*t*lp);
                if(pm > float(Xk)/float(x)) break; pm*=float(x);
            }
        }
    }
    return N0f(t) - acc/PI;
}

// locate zero #k (1-based) by bisection; bracket seeded from the smooth
// estimate (independent => all new zeros locate in parallel this frame)
float locateSmooth(int k, int Xknown, float kmin){
    float tgt=float(k)-0.5, g=smoothGamma(k);
    float Xkf=pow(g/TWO_PI, kmin/TWO_PI);                    // per-zero GHY cutoff
    int Xk=int(min(Xkf, float(Xknown))); if(Xk<2)Xk=2;
    float gap=TWO_PI/log(g/TWO_PI);
    for(float w=1.5; w<=9.0; w+=1.5){
        float lo=g-w*gap, hi=g+w*gap;
        float flo=FB(lo,Xk)-tgt, fhi=FB(hi,Xk)-tgt;
        if(flo*fhi<=0.0){
            for(int it=0;it<50;it++){ float mid=0.5*(lo+hi), fm=FB(mid,Xk)-tgt;
                if(flo*fm<=0.0){hi=mid;fhi=fm;} else {lo=mid;flo=fm;} }
            return 0.5*(lo+hi);
        }
    }
    return g;
}

void mainImage(out vec4 O, in vec2 F){
    ivec2 P=ivec2(F);
    int i = P.x + P.y*STRIDE;

    if(iFrame==0){                                          // seed the state
        float zr=0.0, gi=0.0, b=0.0, a=0.0;
        if(i<NSEED) zr=seed(i);
        if(i==0){ b=float(NSEED); a=1.0; }                  // nz=40, Xknown=1
        else if(i==1){ b=KFIX; a=0.0; }                     // kmin=5, iter=0
        O=vec4(zr,gi,b,a); return;
    }

    vec4 c0=texelFetch(iChannel0, idx2(0),0);
    vec4 c1=texelFetch(iChannel0, idx2(1),0);
    int   nz=int(c0.b+0.5), Xknown=int(c0.a+0.5);
    float kmin=c1.b; int iter=int(c1.a+0.5);

    // ---- plan this frame (deterministic; every texel computes the same plan) ----
    int reach=int(float(nz)/2.7);
    int new_X=Xknown+BAND; if(new_X>reach)new_X=reach; if(new_X<Xknown)new_X=Xknown;
    int new_nz=nz;
    if(new_X>=8){
        int tgt=nz+BATCH, nfr=nFrontier(new_X,kmin);
        if(tgt>nfr)tgt=nfr; if(tgt>NLIM)tgt=NLIM;
        if(tgt>nz) new_nz=tgt;          // only ever GROW the zero list (never regress)
    }

    vec4 old=texelFetch(iChannel0, P, 0);
    float zr=old.r, gi=old.g, b=old.b, a=old.a;

    if(i==0){ b=float(new_nz); a=float(new_X); }            // control 0
    else if(i==1){ b=kmin; a=float(iter+1); }               // control 1

    if(i<NLIM && i>=nz && i<new_nz)                         // newly located zero
        zr = locateSmooth(i+1, Xknown, kmin);

    if(i>=2 && i<=new_X && i>Xknown){                       // newly resolved integer
        float r=detFire(float(i), nz);
        gi = (r>DTAU && !isPP(i)) ? 1.0 : 0.0;
    }

    O=vec4(zr,gi,b,a);
}
