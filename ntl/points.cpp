#include <NTL/tools.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2E.h>

NTL_CLIENT

uint n;
uint N;
uint s;
GF2X fbar;
ZZ* pow2;
ZZX plift;
ZZX* pmod;
GF2E* preRoots;
ZZX* preFrobenius;
ZZ* preTraces;
ZZ exp4;
ZZX tmp[1];

GF2X to_GF2X(const ZZX& a, uint n) {
  int i;
  GF2X b;
  for (i = n - 1; i >= 0; i--) {
    SetCoeff(b, i, IsOdd(coeff(a, i)));
  }
  return b;
}

GF2E to_GF2E(const ZZX& a, uint n) {
  return to_GF2E(to_GF2X(a, n));
}

ZZX to_ZZX(const GF2X& a, uint n) {
  int i;
  ZZX x;
  for (i = n - 1; i >= 0; i--) {
    SetCoeff(x, i, IsOne(coeff(a, i)));
  }
  return x;
}

ZZX to_ZZX(const GF2E& a,  uint n) {
  return to_ZZX(rep(a), n);
}

ZZX operator%(const ZZX& x, const ZZ& a) {
  int i;
  uint n = deg(x) + 1;
  ZZX y(INIT_SIZE, n);
  for (i = n - 1; i >= 0; i-- ) {
    SetCoeff(y, i, coeff(x, i) % a);
  }
  return y;
}

ZZX Mod(const ZZX& x, const ZZ& a, uint n) {
  int i;
  ZZX y(INIT_SIZE, n);
  for (i = n - 1; i >= 0; i-- ) {
    SetCoeff(y, i, coeff(x, i) % a);
  }
  return y;
}

ZZX Unitary(const ZZX& f, uint n, uint N) {
  if (IsOne(LeadCoeff(f))) {
    return ZZX(f);
  }
  return Mod(-f, to_ZZ(1) << N, n + 1);
}

GF2E pthRoot(const GF2E& a) {
  int j, k, l;
  GF2E b;
  GF2X c(INIT_SIZE, n);
  for (j = 0; j < 2; j++) {
    c = 0;
    k = 0;
    l = j;
    while (l < n) {
      SetCoeff(c, k, IsOne(coeff(rep(a), l)));
      k++;
      l += 2;
    }
    b += to_GF2E(c)*(preRoots[j]);
  }
  return b;
}

ZZX qadicInv(const ZZX& x, uint M) {
  if (M == 1) {
    return to_ZZX(1/to_GF2E(x, n), n);
  }
  ZZX& pmodM = (M > N)? pmod[N] : pmod[M];
  uint O = (M >> 1) + (M % 2);
  ZZX y = qadicInv(x, O);
  ZZX xy = MulMod(Mod(x, pow2[M], n), y, pmodM);
  return Mod(MulMod(y, 2 - xy, pmodM), pow2[M], n);
}

ZZ qadicTrace(const ZZX& x, uint N) {
  int i;
  ZZ t = to_ZZ(0);
  for (i = 0; i < n; i++) {
    t += coeff(x, i)*preTraces[i];
  }
  return t % pow2[N];
}

uint valuation(const ZZX& x, uint N) {
  int i;
  uint v = N;
  uint u;
  ZZ c;
  for (i = 0; i < n; i++) {
    u = 0;
    c = ZZ(coeff(x, i));
    while (c % 2 != 1) {
      u++;
      c >>= 1;
    }
    if (u < v) {
      v = u;
    }
  }
  return v;
}

ZZX Frobenius(const ZZX& x, uint N) {
  int i;
  ZZX y;
  for (i = 0; i < n; i++) {
    y += coeff(x, i)*Mod(preFrobenius[i], pow2[N], n);
  }
  return Mod(y, pow2[N], n);
}

ZZX TeichmullerFieldPolInc2(const ZZX& g0, const ZZX& g1, const ZZX& V, uint N) {
  int i;
  if (N == 1) {
    return Mod(V, pow2[1], n + 1);
  }
  uint M = (N >> 1) + (N % 2);
  ZZX e = TeichmullerFieldPolInc2(g0, g1, V, M);
  ZZX e0 = ZZX(INIT_SIZE, (n >> 1) + 1);
  ZZX e1 = ZZX(INIT_SIZE, (n + 1) >> 1);
  for (i = n; i >= 0; i--) {
    if (i % 2) {
      // e1(x^2) = (e(x) - e(-x)) / (2 x)
      SetCoeff(e1, i >> 1, coeff(e, i));
    } else {
      // e0(x^2) = (e(x) + e(-x)) / 2
      SetCoeff(e0, i >> 1, coeff(e, i));
    }
  }
  // U = (V + e - 2(g0 e0 - x g1 e1)) / 2^M
  ZZX U = Mod(V + e - 2*(g0*e0 - ((g1*e1) << 1)), pow2[N], n + 1) / pow2[M];
  ZZX D = TeichmullerFieldPolInc2(g0, g1, U, N - M);
  //cout << "e = " << e << endl;
  //cout << "V = " << V << endl;
  //cout << "U = " << U << endl;
  //cout << "D = " << D << endl;
  return Mod(e + pow2[M]*D, pow2[N], n + 1);
}

ZZX TeichmullerFieldPol2(const ZZX& f, uint N) {
  int i;
  if (N == 1) {
    return ZZX(f);
  }
  uint M = (N >> 1) + (N % 2);
  ZZX g = TeichmullerFieldPol2(f, M);
  //cout << "g[" << M << "] = " << g << endl;
  ZZX g0 = ZZX(INIT_SIZE, (n >> 1) + 1);
  ZZX g1 = ZZX(INIT_SIZE, (n + 1) >> 1);
  for (i = n; i >= 0; i--) {
    if (i % 2) {
      // g1(x^2) = (g(x) - g(-x)) / (2 x)
      SetCoeff(g1, i >> 1, coeff(g, i));
    } else {
      // g0(x^2) = (g(x) + g(-x)) / 2
      SetCoeff(g0, i >> 1, coeff(g, i));
    }
  }
  // V = (g - g0^2 - x g1^2) / 2^M
  ZZX V = Mod(g - sqr(g0) + (sqr(g1) << 1), pow2[N], n + 1) / pow2[M];
  //cout << "V[" << M << "] = " << V << endl;
  ZZX d = TeichmullerFieldPolInc2(g0, g1, V, N - M);
  //cout << "d[" << M << "] = " << d << endl;
  return Mod(g + pow2[M]*d, pow2[N], n + 1);
}

ZZX EvalG2(const ZZX& x, const ZZX& y, uint N) {
  // (X+2Y+8XY)^2 + Y + 4XY
  ZZX xy = MulMod(x, y, pmod[N]);
  return Mod(SqrMod(x + 2*y + 8*xy, pmod[N]) + y + 4*xy, pow2[N], n);
}

ZZX EvalG2x(const ZZX& x, const ZZX& y, uint N) {
  // 2(X+2Y+8XY)(1+8Y) + 4Y
  return Mod(2*MulMod(x + 2*y + 8*MulMod(x, y, pmod[N]), 1 + 8*y, pmod[N]) + 4*y, pow2[N], n);
}

ZZX EvalG2y(const ZZX& x, const ZZX& y, uint N) {
  // (4(X+2Y+8XY)+1)(1+4X)
  return Mod(MulMod(4*(x + 2*y + 8*MulMod(x, y, pmod[N])) + 1, 1 + 4*x, pmod[N]), pow2[N], n);
}

ZZX ArtinSchreierRootII(const ZZX& a, const ZZX& b, const ZZX& c, uint N) {
  if (N == 1) {
    //cout << "inv = " << 1/to_GF2E(a, n) << endl;
    //cout << "mul = " << to_GF2E(c, n)/to_GF2E(a, n) << endl;
    //cout << "neg = " << -to_GF2E(c, n)/to_GF2E(a, n) << endl;
    //cout << "root = " << to_ZZX(pthRoot(-to_GF2E(c, n)/to_GF2E(a, n)), n) << endl;
    return to_ZZX(pthRoot(-to_GF2E(c, n)/to_GF2E(a, n)), n);
  }
  uint M = (N >> 1) + (N % 2);
  uint O = N - M;
  ZZX y = ArtinSchreierRootII(a, b, c, M);
  //cout << "y = " << y << endl;
  ZZX d = Mod((a*Frobenius(y, N) + b*y + c) % pmod[N], pow2[N], n) / pow2[M];
  //cout << "d = " << d << endl;
  ZZX E = ArtinSchreierRootII(a, b, d, O);
  //cout << "E = " << E << endl;
  return Mod(y + pow2[M]*E, pow2[N], n);
}

ZZX GenNewtonLiftII(const ZZX& x0, uint N) {
  if (N == 1) {
    return ZZX(x0);
  }
  uint M = (N >> 1) + (N % 2);
  ZZX x = GenNewtonLiftII(x0, M);
  //cout << "N = " << N << endl;
  //cout << "M = " << M << endl;
  //cout << "x = " << x << endl;
  ZZX y = Frobenius(x, N);
  //cout << "y = " << y << endl;
  ZZX V = EvalG2(x, y, N) / pow2[M];
  //cout << "V = " << V << endl;
  ZZX Dx = EvalG2x(x, y, M);
  //cout << "Dx = " << Dx << endl;
  ZZX Dy = EvalG2y(x, y, M);
  //cout << "Dy = " << Dy << endl;
  ZZX E = ArtinSchreierRootII(Dy, Dx, V, M);
  //cout << "E = " << Mod(pow2[M]*E, pow2[N], n) << endl;
  return Mod(x + pow2[M]*E, pow2[N], n);
}

ZZ NormChar2(const ZZX& a, uint N, uint s) {
  int i, j;
  ZZX z = ZZX(a);
  //cout << "s = " << s << endl;
  for (i = 0; i < s; i++) {
    z = Mod(SqrMod(z, pmod[N]), pow2[N + s], n);
  }
  z -= 1;
  //cout << "a = " << a << endl;
  z /= 2;
  //cout << "z = " << z << endl;
  uint v = valuation(z, N + s - 1);
  ZZX c = Mod(MulMod(z, qadicInv(1 + z, N + s - 1), pmod[N]), pow2[N + s - 1], n);
  ZZX c2 = Mod(SqrMod(c, pmod[N]), pow2[N + s - 1], n);
  //cout << "c = " << c << endl;
  //cout << "c2 = " << c2 << endl;

  ZZX log_a;
  i = 1;
  j = v;
  while (j < N + s - 1) {
    tmp[0] = c * InvMod(to_ZZ(i), pow2[N + s - 1]);
    log_a += tmp[0];
    c = Mod(MulMod(c, c2, pmod[N]), pow2[N + s - 1], n);
    i += 2;
    j += 2*v;
  }
  log_a = Mod(log_a, pow2[N + s - 1], n);

  ZZX w = log_a / pow2[s - 1];
  //cout << "w = " << w << endl;

  ZZ u = qadicTrace(w, N) / 4;
  //cout << "u = " << u << endl;

  return PowerMod(exp4, u, pow2[N]);
}

void GF2nInit() {
  fbar = BuildSparseIrred_GF2X(n);
  cout << "fbar = " << fbar << endl;
  GF2E::init(fbar);
}

void PrecompPow2(uint N, uint s) {
  int i;
  pow2 = new ZZ[N + s + 1];
  pow2[0] = 1;
  for (i = 1; i <= N + s; i++) {
    pow2[i] = pow2[i - 1] << 1;
  }
}

void PrecompMod(uint N) {
  int i;
  pmod = new ZZX[N + 1];
  for (i = 0; i <= N; i++) {
    pmod[i] = Mod(plift, pow2[i], n + 1);
  }
}

void PrecompRoots(uint N) {
  int i;
  GF2X x;
  SetX(x);
  GF2E t = power(to_GF2E(x), to_ZZ(1) << (n - 1));
  preRoots = new GF2E[2];
  preRoots[0] = 1;
  //cout << "Roots[" << 0 << "] = " << preRoots[0] << endl;
  for (i = 1; i < 2; i++) {
    preRoots[i] = power(t, i);
    //cout << "Roots[" << i << "] = " << preRoots[i] << endl;
  }
}

void PrecompFrobenius(uint N) {
  int i;
  ZZX x2 = to_ZZX(1) << 2;
  preFrobenius =  new ZZX[n];
  preFrobenius[0] = 1;
  //cout << "F[0] = " << preFrobenius[0] << endl;
  for (i = 1; i < n; i++) {
    preFrobenius[i] = Mod(MulMod(preFrobenius[i - 1], x2, pmod[N]), pow2[N], n);
    //cout << "F[" << i << "] = " << preFrobenius[i] << endl;
  }
}

void PrecompTraces(uint N) {
  int i, j;
  preTraces = new ZZ[n];
  preTraces[0] = to_ZZ(n);
  ZZ t;
  for (i = 1; i < n; i++) {
    t = -i*coeff(plift, n-i);
    for (j = 1; j < i; j++) {
      t -= coeff(plift, n-j)*preTraces[i-j];
    }
    preTraces[i] = t % pow2[N];
  }
}

void PrecompExp4(uint N) {
  int i = 2;
  exp4 = to_ZZ(1);
  ZZ x = to_ZZ(1);
  ZZ y = to_ZZ(1);
  for (i = 1; i < N; i++) {
    x <<= 2;
    y = MulMod(y, i, pow2[N]);
    x >>= MakeOdd(y);
    exp4 = AddMod(exp4, MulMod(x, InvMod(y, pow2[N]), pow2[N]), pow2[N]);
  }
}

void PrecompQAdics(uint N) {
  PrecompMod(N);
  PrecompRoots(N);
  PrecompFrobenius(N);
  PrecompTraces(N);
  PrecompExp4(N);
}

int main(int argc, char* argv[]) {
  double time;

  // Extension degree
  n = atoi(argv[1]);
  GF2nInit();

  // Precision
  N = (n >> 1) + (n % 2) + 2;
  // Extra precision for norm
  s = to_uint(to_ZZ(sqrt(to_RR(N)))) >> 1;
  // Precompute powers of 2
  PrecompPow2(N, s);
  cout << "Needed precision " << N << endl;

  // Teichmüller lift
  time = GetTime();
  plift = Unitary(TeichmullerFieldPol2(to_ZZX(fbar, n + 1), N), n + 1, N);
  cout << "Teichmuller lift done in ";
  PrintTime(cout, GetTime() - time);
  cout << " seconds." << endl;
  cout << "plift = " << plift << endl;

  // Precomputations for q-adics
  time = GetTime();
  PrecompQAdics(N);
  cout << "Precomputations done in ";
  PrintTime(cout, GetTime() - time);
  cout << " seconds." << endl;

  // random elliptic curve
  ZZX a = ZZX(n - 2, 1);
  cout << "a = " << a << endl;

  time = GetTime();
  ZZX x = GenNewtonLiftII(a, N);
  cout << "Lift done in ";
  PrintTime(cout, GetTime() - time);
  cout << " seconds." << endl;
  cout << "x = " << x << endl;

  time = GetTime();
  cout << "t = " << qadicInv(1 + 4*x, N) << endl;
  ZZ t = NormChar2(qadicInv(1 + 4*x, N), N, s);
  cout << "Norm done in ";
  PrintTime(cout, GetTime() - time);
  cout << " seconds." << endl;
  cout << "t = " << t << endl;

  if (sqr(t) > (to_ZZ(1) << (n + 2))) {
    t -= to_ZZ(1) << N;
  }
  cout << "t = " << t << endl;

  ZZ c = (to_ZZ(1) << n) + 1 - t;
  cout << "#E_a = " << c << endl;

  return 0;
}
