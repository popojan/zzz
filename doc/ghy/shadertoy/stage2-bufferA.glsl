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

// detector r(x) from the nz committed zeros (.r channel of Buffer A): the reading at six
// truncations of the zero list (0.5 .. 1.0 of nz; prefix sums of one pass), each scaled by
// its own top zero, and the median of the six.  A single reading oscillates with the
// truncation height, so near the reach a composite next to a prime can ride a sidelobe.
float detFire(float x, int nz){
    float lx=log(x), s=0.0, base=1.0 + 1.0/(x - x*x*x), c4=4.0/sqrt(x), v[6];
    int b=0, nb=nz/2;
    for(int k=0;k<NLIM;k++){
        if(k>=nz) break;
        float zk=texelFetch(iChannel0, idx2(k),0).r;
        s += pcos(zk*lx);
        if(k+1==nb){
            v[b]=(base - c4*s)/peakf(x, zk);
            b++; nb = (b==5) ? nz : int(floor((0.5+0.1*float(b))*float(nz)));
        }
    }
    for(int i=1;i<6;i++){ float t=v[i]; int j=i-1;            // insertion sort
        for(int q=0;q<6;q++){ if(j<0 || v[j]<=t) break; v[j+1]=v[j]; j--; }
        v[j+1]=t; }
    return 0.5*(v[2]+v[3]);
}

// F_B(t): method-B count from the prime-indicator field (.g), primes <= Xk
float FB(float t, int Xk){
    float acc=0.0, w=BOXC*TWO_PI/log(t/TWO_PI);          // box average over [t-w, t+w]
    for(int x=2;x<XKMAX;x++){
        if(x>Xk) break;
        if(texelFetch(iChannel0, idx2(x),0).g > 0.5){       // x is prime
            float lp=log(float(x)), pm=float(x);
            for(int m=1;m<20;m++){
                float hf=pow(float(x), -0.5*float(m));
                float a=float(m)*w*lp;                      // sinc weight: damps the truncation edge
                acc += (hf/float(m))*psin(float(m)*t*lp)*(a>0.0 ? sin(a)/a : 1.0);
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
    // decide x only where the LOWEST truncation still resolves unit spacing; below about
    // gamma(nz/2) = 1.8*x a composite between twin primes reads as a prime
    int reach=int(texelFetch(iChannel0, idx2(nz/2-1),0).r/UMIN);
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
