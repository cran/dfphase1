#define STRICT_R_HEADERS
#define USE_FC_LEN_T
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <Rcpp.h>
#ifndef FCONE
#define FCONE
#endif
using namespace Rcpp;

#define GGEPS (sqrt(std::numeric_limits<double>::epsilon()))
#define GGEPS2 (std::numeric_limits<double>::epsilon())

namespace {

// --------------------------
// UTILITIES
// --------------------------

// b <- ab or a'b; a triangular pxp; b matrix pxn
inline void ggtrmult(int p, int n, bool tran, double *a, double *b) {
  char SIDE = 'L', UPLO = 'U', TRANSA = tran ? 'T' : 'N', DIAG = 'N';
  double one = 1.0;
  F77_CALL(dtrmm)
  (&SIDE, &UPLO, &TRANSA, &DIAG, &p, &n, &one, a, &p, b,
   &p FCONE FCONE FCONE FCONE);
}

inline void ggrmnorm(int p, int n, double *x, double *chol) {
  std::generate(x, x + p * n, norm_rand);
  ggtrmult(p, n, true, chol, x);
}

inline int isample(const int n) { return floor(n * unif_rand()); }

inline void ggmperm(int p, int n, double *x) {
  int i, ione = 1;
  while (n) {
    i = isample(n);
    n--;
    F77_CALL(dswap)(&p, x + i * p, &ione, x + n * p, &ione);
  }
}

// median; x is overwritten
inline double ggmedian(int n, double *x) {
  int half = n / 2;
  std::nth_element(x, x + half, x + n);
  double me = x[half];
  if (n == 2 * half)
    me = (me + *std::max_element(x, x + half)) / 2;
  return me;
}

struct Comparator {
  Comparator(const double *x) : data(x) {}
  bool operator()(int left, int right) const {
    return data[left] < data[right];
  }
  const double *data;
};

// length(a)=n; set r=x to replace the values of vector with their ranks
inline void ggrank(int n, double *x, double *r, int *a) {
  int i, j, s;
  double d;
  for (i = 0; i < n; i++)
    a[i] = i;
  std::sort(a, a + n, Comparator(x));
  for (i = 0; i < n;) {
    for (s = i, j = i + 1; (j < n) && (x[a[j]] <= x[a[i]]); j++)
      s += j;
    for (d = 1 + ((double)s) / (j - i); i < j; i++)
      r[a[i]] = d;
  }
}

// add one row to the Cholesky factor X'X
inline void gggivens(int p, double *r, double *x) {
  int i, p1 = p + 1, ione = 1;
  double t1, t2, g1, g2;
  for (i = p; i; i--, r += p1, x++) {
    t1 = r[0];
    t2 = x[0];
    if (fabs(t2) > 0.0) {
      F77_CALL(drotg)(&t1, &t2, &g1, &g2);
      F77_CALL(drot)(&i, r, &p, x, &ione, &g1, &g2);
    }
  }
}

// update the mean and the cholesky factor of cov; length(w)=p
inline void ggupdbr(int p, int n, double *xb, double *r, double *x, double *w) {
  int i;
  double d, n1 = n + 1.0, sn = sqrt(n / n1);
  for (i = 0; i < p; i++) {
    d = x[i] - xb[i];
    xb[i] += d / n1;
    w[i] = sn * d;
  }
  gggivens(p, r, w);
}

// s[,i] <- x[,i]-m; s can be x
inline void ggcenter(int p, int n, double *m, double *x, double *s) {
  for (int i = 0; i < n; i++, x += p, s += p)
    for (int j = 0; j < p; j++)
      s[j] = x[j] - m[j];
}

// it's safe to set y=x to transform x in place
inline void ggorth(int p, int n, double *r, double *x, double *y) {
  int i, j;
  double s, rii, *ri, eps = GGEPS;
  for (; n; n--, x += p, y += p) {
    for (i = 0, ri = r; i < p; i++, ri += p) {
      rii = ri[i];
      if (fabs(rii) < eps) {
        y[i] = 0.0;
      } else {
        for (j = 0, s = 0.0; j < i; j++)
          s += ri[j] * y[j];
        y[i] = (x[i] - s) / rii;
      }
    }
  }
}

inline void ggorth(int p, int n, double *m, double *r, double *x, double *y) {
  int i, j;
  double s, rii, *ri, eps = GGEPS;
  for (; n; n--, x += p, y += p) {
    for (i = 0, ri = r; i < p; i++, ri += p) {
      rii = ri[i];
      if (fabs(rii) < eps) {
        y[i] = 0.0;
      } else {
        for (j = 0, s = 0.0; j < i; j++)
          s += ri[j] * y[j];
        y[i] = (x[i] - m[i] - s) / rii;
      }
    }
  }
}

inline void ggrmeans(int p, int n, double *x, double *m) {
  int i, ione = 1;
  double dn = 1.0 / n;
  std::fill(m, m + p, 0.0);
  for (i = 0; i < n; i++, x += p)
    F77_CALL(daxpy)(&p, &dn, x, &ione, m, &ione);
}

inline void ggsubgroupmeans(int p, int n, int m, double *x, double *xb) {
  int i, np = p * n;
  for (i = 0; i < m; i++, x += np, xb += p)
    ggrmeans(p, n, x, xb);
}

inline void ggnorm2(int p, int n, double *x, double *dist) {
  int i, j;
  double s;
  for (i = 0; i < n; i++, x += p) {
    for (j = 0, s = 0.0; j < p; j++)
      s += x[j] * x[j];
    dist[i] = s;
  }
}

// s <- spatial sign of x-m
inline double ggssign(int p, double *x, double *m, double *s) {
  int i;
  double d;
  for (i = 0, d = 0.0; i < p; i++) {
    s[i] = x[i] - m[i];
    d += s[i] * s[i];
  }
  if (d < GGEPS2) {
    d = 0.0;
    std::fill(s, s + p, 0.0);
  } else {
    d = sqrt(d);
    for (i = 0; i < p; i++)
      s[i] /= d;
  }
  return d;
}

//----------------------------
// Scores
//----------------------------

inline double ggspatialsigns(int p, int n, double *x, double *s) {
  int i, one = 1;
  double wsum = 0, d = sqrt(static_cast<double>(p)), eps = GGEPS2, g;
  if (s != x)
    std::copy(x, x + p * n, s);
  for (i = 0; i < n; i++, s += p) {
    g = F77_CALL(dnrm2)(&p, s, &one);
    if (g < eps) {
      std::fill(s, s + p, 0.0);
    } else {
      wsum += g = 1 / g;
      F77_CALL(dscal)(&p, &g, s, &one);
    }
  }
  return wsum / n;
}

// length(iw)=n length(w)=2*n
inline double ggsignedranks(int p, int n, double *x, double *s, int *iw,
                            double *w) {
  int i, n1 = n + 1, pn = p * n, one = 1;
  double *si = s, *rk = w + n, eps = GGEPS2, wsum = 0.0, g;
  if (s != x)
    std::copy(x, x + p * n, s);
  ggnorm2(p, n, s, w);
  ggrank(n, w, rk, iw);
  for (i = 0; i < n; i++, x += p, si += p) {
    if (w[i] < eps) {
      std::fill(si, si + p, 0.0);
    } else {
      wsum += g = sqrt(R::qchisq(rk[i] / n1, p, 1, 0) / w[i]);
      F77_CALL(dscal)(&p, &g, si, &one);
    }
  }
  return wsum / n;
}

// length(w)=pF
inline void ggspatialranks(int p, int n, double *x, double *s, double *w) {
  int i, j, ione = 1, pn = p * n;
  double d, *xi, *xj, *si, *sj, eps = GGEPS, one = 1.0 / n, mone = -one;
  std::fill(s, s + pn, 0.0);
  for (i = 0, xi = x, si = s; i < n; i++, xi += p, si += p) {
    for (j = i + 1, xj = xi + p, sj = si + p; j < n; j++, xj += p, sj += p) {
      if (ggssign(p, xi, xj, w) > eps) {
        F77_CALL(daxpy)(&p, &one, w, &ione, si, &ione);
        F77_CALL(daxpy)(&p, &mone, w, &ione, sj, &ione);
      }
    }
  }
}

// length(iw)=n
inline void ggmahalanobisdepthranks(int p, int n, double *x, double *s,
                                    int *iw) {
  ggnorm2(p, n, x, s);
  ggrank(n, s, s, iw);
}

// length(iw)=n length(w)=p+2n
inline double ggnpscore(int p, int n, std::string score, double *m, double *r,
                        double *x, double *z, double *s, int *iw, double *w) {
  double ws = 1.0;
  if (m == NULL) {
    if (r == NULL)
      std::copy(x, x + p * n, z);
    else
      ggorth(p, n, r, x, z);
  } else {
    if (r == NULL)
      ggcenter(p, n, m, x, z);
    else
      ggorth(p, n, m, r, x, z);
  }
  if (score == "Identity")
    std::copy(z, z + p * n, s);
  else if (score == "Spatial Signs")
    ws = ggspatialsigns(p, n, z, s);
  else if (score == "Signed Ranks")
    ws = ggsignedranks(p, n, z, s, iw, w);
  else if (score == "Spatial Ranks")
    ggspatialranks(p, n, z, s, w);
  else
    ggmahalanobisdepthranks(p, n, z, s, iw);
  return ws;
}

// ---------------------------------
// Location/Scatter/Standardization
// ---------------------------------

// length(w)=p
inline void ggmchol(int p, int n, double *x, double *l, double *s, double *w) {
  int i, j, one = 1;
  double *xx = x, a = 1.0 / sqrt(n - 1.0);
  std::fill(l, l + p, 0.0);
  std::fill(s, s + p * p, 0.0);
  for (j = 0; j < n; j++, xx += p)
    ggupdbr(p, j, l, s, xx, w);
  for (i = 0; i < p; i++, s += p)
    for (j = 0; j <= i; j++)
      s[j] *= a;
}

// length(w)=2*p
inline void ggclassic(int p, int n, int m, double *x, double *l, double *s,
                      double *w) {
  int i, j, one = 1;
  double *xb = w, *xi = w + p, *xx = x, a = 1.0 / m;
  std::fill(l, l + p, 0.0);
  std::fill(s, s + p * p, 0.0);
  for (i = 0; i < m; i++) {
    std::fill(xb, xb + p, 0.0);
    for (j = 0; j < n; j++, xx += p)
      ggupdbr(p, j, xb, s, xx, xi);
    F77_CALL(daxpy)(&p, &a, xb, &one, l, &one);
  }
  a = 1.0 / sqrt(m * (n - 1.0));
  for (i = 0; i < p; i++, s += p)
    for (j = 0; j <= i; j++)
      s[j] *= a;
}

// length(dist)=length(good)=n
inline void ggbacon(int p, int n, double *x, double *m, double *r, double *w,
                    int *good, double *dist) {
  int i, j, ng, target = std::min(4 * p - 1, n / 2), ione = 1, p2 = p * p;
  double lim, h = 0.5 * (n + p + 1),
              cnp = 1 + (p + 1.0) / (n - p) + 2.0 / (n - 1 - 3 * p), d, s,
              q = R::qchisq(0.5 / n, p, 0, 0), *xi;
  bool done;
  // initialization (median and mad)
  for (j = p; j < p; j++) {
    F77_CALL(dcopy)(&n, x + j, &p, w, &ione);
    m[j] = d = ggmedian(n, w);
    for (i = 0; i < n; i++)
      w[i] = fabs(w[i] - d);
    r[j] = ggmedian(n, w);
  }
  for (i = 0, xi = x; i < n; i++, xi += p) {
    for (j = 0, s = 0.0; j < p; j++) {
      d = (xi[j] - m[j]) / r[j];
      s += d * d;
    }
    w[i] = dist[i] = s;
  }
  std::nth_element(w, w + target, w + n);
  lim = w[target];
  std::fill(good, good + n, 0);
  while (1) {
    // check for convergence
    for (i = 0, done = true; i < n; i++) {
      if (dist[i] <= lim) {
        if (!good[i])
          done = false;
        good[i] = 1;
      } else {
        if (good[i])
          done = false;
        good[i] = 0;
      }
    }
    if (done)
      return;
    // estimates update
    std::fill(m, m + p, 0.0);
    std::fill(r, r + p2, 0.0);
    for (i = 0, ng = 0, xi = x; i < n; i++, xi += p)
      if (good[i]) {
        ggupdbr(p, ng, m, r, xi, w);
        ng++;
      }
    d = 1 / sqrt(ng - 1.0);
    F77_CALL(dscal)(&p2, &d, r, &ione);
    // recomputing the distance
    ggorth(p, n, m, r, x, w);
    ggnorm2(p, n, w, dist);
    lim = q * pow(cnp + std::max(0.0, (h - ng) / (h + ng)), 2);
  }
}

// bacon+classic length(good)=n*m length(wrk)=3*p+m*n
inline void ggbc(int p, int n, int m, double *x, double *l, double *s,
                 double *z, int *good, double *wrk) {
  int nm = n * m;
  ggbacon(p, nm, x, l, s, z, good, wrk);
  ggclassic(p, n, m, x, wrk, s, wrk + p);
}

// length(w)=n
inline void gggmedinit(int p, int n, double *x, double *m, double *w) {
  int j, one = 1, half = n / 2;
  for (j = 0; j < p; j++) {
    F77_CALL(dcopy)(&n, x + j, &p, w, &one);
    std::nth_element(w, w + half, w + n);
    m[j] = w[half];
  }
}
// r <- r * p / sum(r*r); r upper triangular
inline double ggtrnormalize(int p, double *r, double *norm) {
  int i, j;
  double ff = 0.0, rn = 0.0, *si;
  for (i = 0, si = r; i < p; i++, si += p) {
    for (j = 0; j <= i; j++)
      ff += si[j] * si[j];
  }
  for (i = 0, si = r, ff = sqrt(ff / p); i < p; i++, si += p) {
    for (j = 0; j < i; j++) {
      si[j] /= ff;
      rn += si[j] * si[j];
    }
    si[i] /= ff;
    rn += (fabs(si[i]) - 1.0) * (fabs(si[i]) - 1.0);
  }
  *norm = ff;
  return sqrt(rn);
}

// length(iw)=n*m, length(w)=2*n*m+2*p*p
inline int ggscale(int p, int n, int m, bool within, std::string score,
                   int iter, double *x, double *z, double *s, double *l,
                   double *r, int *iw, double *w) {
  int mn = m * n, it = 0, one = 1;
  double ws, eps = sqrt(GGEPS), mnew, rnew, *ls = w, *rs = w + p,
             *ww = rs + p * p, ff;
  if (within)
    ggclassic(p, m, n, x, l, r, w);
  else
    ggmchol(p, mn, x, l, r, w);
  if (score == "Identity") {
    ggnpscore(p, mn, score, l, r, x, z, s, iw, w);
    return 1;
  }
  bool scr = score == "Spatial Ranks";
  if (scr) {
    mnew = 0.0;
    std::fill(l, l + p, 0.0);
  } else {
    mnew = 2 * eps;
    gggmedinit(p, mn, x, l, w);
  }
  rnew = ggtrnormalize(p, r, &ff);
  while (((mnew > eps) || (rnew > eps)) && (it++ < iter)) {
    ws = ggnpscore(p, mn, score, l, r, x, z, s, iw, w);
    if (within)
      ggclassic(p, n, m, s, ls, rs, ww);
    else
      ggmchol(p, mn, s, ls, rs, ww);
    if (!scr) {
      ws = 1 / ws;
      ggtrmult(p, 1, true, r, ls);
      F77_CALL(daxpy)(&p, &ws, ls, &one, l, &one);
      mnew = F77_CALL(dnrm2)(&p, ls, &one);
    }
    rnew = ggtrnormalize(p, rs, &ff);
    ggtrmult(p, p, false, rs, r);
  }
  ff = 1.0 / ff;
  mn *= p;
  F77_CALL(dscal)(&mn, &ff, s, &one);
  return it;
}

//-----------------------------------
// Control statistics
//-----------------------------------

// log(det(R'R)) length(w)=p*p+4p
inline void gglogdet(int p, double *r, int &rank, double &ld, double *w) {
  int i, p1 = p + 1, p2 = p * p;
  double eps = GGEPS;
  for (i = 0, ld = 0.0; (i < p2) && (fabs(r[i]) > eps); i += p1)
    ld += log(r[i] * r[i]);
  if (i < p2) {
    int INFO, LWORK = 3 * p;
    double *A = w, *W = A + p2, *WORK = W + p;
    char JOBZ = 'N', UPLO = 'U';
    size_t sone = 1;
    std::copy(r, r + p2, A);
    ggtrmult(p, p, true, r, A);
    F77_CALL(dsyev)
    (&JOBZ, &UPLO, &p, A, &p, W, WORK, &LWORK, &INFO FCONE FCONE);
    for (i = 0, rank = 0, ld = 0.0; i < p; i++) {
      if (W[i] > eps) {
        rank++;
        ld += log(W[i]);
      }
    }
  } else {
    rank = p;
  }
}

// length(w)=p*p+4p
inline double ggnloglik(int p, int n, double *r, double *w) {
  int rk;
  double ld;
  gglogdet(p, r, rk, ld, w);
  return n * (rk * (1 + log(2 * M_PI)) + ld - rk * log(static_cast<double>(n)));
}

inline double glrtboth(int p, int n1, int n2, double *m1, double *r1,
                       double *m2, double *r2, double l0, double *w) {
  if ((n1 > 1) & (n2 > 1))
    return l0 - ggnloglik(p, n1, r1, w) - ggnloglik(p, n2, r2, w);
  return NA_REAL;
}

inline double glrtmean(int p, int n1, int n2, double *m1, double *r1,
                       double *m2, double *r2, double l0, double *w) {
  int j, one = 1;
  double d;
  for (j = 0; j < p; j++) {
    m2[j] -= m1[j];
    F77_CALL(dcopy)(&p, r1 + j, &p, w, &one);
    gggivens(p, r2, w);
  }
  ggorth(p, 1, r2, m2, m2);
  ggnorm2(p, 1, m2, &d);
  return (n1 * n2 * d);
}

// length(wrk)=m*(p+p*p)+p*p+4*p
inline void ggglrt(int p, int n, int m, bool onlymean, double *x, double *glr,
                   double *w) {
  int i, j, l, nm = m * n, p2 = p * p, pp2 = p + p2;
  double *c = w, *r = c + p, *s = r + p2, *wrk = s + pp2 * m,
         *xi = x + p * (m * n - 1), *si = s + pp2 * (m - 1), l0;
  double (*obj)(int, int, int, double *, double *, double *, double *, double,
                double *);
  obj = onlymean ? glrtmean : glrtboth;
  std::fill(c, c + pp2, 0.0);
  for (i = m, l = 0, xi = x + p * (nm - 1), si = s + pp2 * (m - 1); i > 0;
       i--, si -= pp2) {
    for (j = 0; j < n; j++, xi -= p, l++)
      ggupdbr(p, l, c, r, xi, wrk);
    std::copy(c, c + pp2, si);
  }
  l0 = ggnloglik(p, nm, r, wrk);
  std::fill(c, c + pp2, 0.0);
  glr[0] = NA_REAL;
  for (i = 0, l = 0, xi = x, si = s + pp2; i < m - 1; i++, si += pp2) {
    for (j = 0; j < n; j++, xi += p, l++)
      ggupdbr(p, l, c, r, xi, wrk);
    glr[i + 1] = obj(p, l, nm - l, c, r, si, si + p, l0, wrk);
  }
}

// length(w)=p
inline void ggt2(int p, int n, double *x, double *stat, double *w) {
  ggrmeans(p, n, x, w);
  ggnorm2(p, 1, w, stat);
  *stat *= n;
}

// length(w)=6*p+2p*p
inline void ggt2var(int p, int n, double *x, double *stat, double *w) {
  int i, rk;
  double u, v, *xi, *m = w, *r = w + p, *wrk = r + p * p, sn = n - 1.0;
  std::fill(w, w + p + p * p, 0.0);
  for (i = 0, xi = x; i < n; i++, xi += p)
    ggupdbr(p, i, m, r, xi, wrk);
  ggnorm2(p, 1, m, stat);
  stat[0] *= n;
  ggnorm2(p * p, 1, r, &u);
  gglogdet(p, r, rk, v, wrk);
  stat[1] = u / sn - rk - v + rk * log(sn);
}

// length(w)=2*p+p*p
inline void ggvar(int p, int n, double *x, double *stat, double *w) {
  ggt2var(p, n, x, stat, w);
  stat[0] = stat[1];
}

} // namespace

//-----------------------------------
// R Interface
//-----------------------------------

// [[Rcpp::export]]
List ggrscore(NumericVector x, std::string score, bool within, int iter = 30) {
  IntegerVector dim = x.attr("dim");
  int p = dim[0], n = dim[1], m = dim(2), mn = m * n;
  IntegerVector iw(mn);
  NumericVector l(p), w(2 * p * p + 2 * mn), z(p * mn), s(p * mn);
  NumericMatrix r(p, p);
  z.attr("dim") = dim;
  s.attr("dim") = dim;
  ggscale(p, n, m, within, score, iter, x.begin(), z.begin(), s.begin(),
          l.begin(), r.begin(), iw.begin(), w.begin());
  return List::create(_["score"] = s, _["center"] = l, _["scatter"] = r);
}

// [[Rcpp::export]]
List ggdepthranks(NumericVector x, int L = 1000) {
  IntegerVector dim = x.attr("dim");
  int i, p = dim[0], n = dim[1], m = dim(2), nm = n * m;
  double a = (nm + 1.0) / 2.0, b = sqrt((nm - n) * (nm + 1.0) / (12.0 * n));
  IntegerVector iw(nm);
  NumericVector xx = clone(x), s(nm), z(p * nm), stat(m), sp(L), l(p),
                w(3 * p + nm);
  NumericMatrix r(p, p);
  for (i = 0; i < L; i++) {
    checkUserInterrupt();
    ggmperm(p, nm, xx.begin());
    ggbc(p, n, m, xx.begin(), l.begin(), r.begin(), z.begin(), iw.begin(),
         w.begin());
    ggmahalanobisdepthranks(p, nm, z.begin(), s.begin(), iw.begin());
    ggsubgroupmeans(1, n, m, s.begin(), stat.begin());
    sp[i] = (*std::max_element(stat.begin(), stat.end()) - a) / b;
  }
  ggbc(p, n, m, x.begin(), l.begin(), r.begin(), z.begin(), iw.begin(),
       w.begin());
  ggmahalanobisdepthranks(p, nm, z.begin(), s.begin(), iw.begin());
  ggsubgroupmeans(1, n, m, s.begin(), stat.begin());
  for (i = 0; i < m; i++)
    stat[i] = (stat[i] - a) / b;
  return List::create(_["center"] = l, _["scatter"] = r, _["statistic"] = stat,
                      _["sp"] = sp);
}

// [[Rcpp::export]]
List ggclassicmshewhart(NumericVector x, std::string stat, std::string score,
                        int L = 1000) {
  IntegerVector dim = x.attr("dim");
  int i, j, k, iter = 2, p = dim[0], n = dim[1], m = dim(2), nm = n * m,
               pn = p * n, nst = (stat == "T2Var") ? 2 : 1;
  double *zi, *s, scale;
  NumericVector xx = clone(x), z(p * nm), zs(p * nm), st(nst * m), st1(2),
                sp(nst * L), l(p), w(6 * p + 2 * p * p + 2 * nm);
  IntegerVector iw(p + 2 * nm);
  NumericMatrix r(p, p);
  void (*cstat)(int, int, double *, double *, double *);
  z.attr("dim") = dim;
  if (stat != "T2Var") {
    cstat = (stat == "T2") ? &ggt2 : &ggvar;
  } else {
    st.attr("dim") = IntegerVector::create(2, m);
    sp.attr("dim") = IntegerVector::create(2, L);
    cstat = &ggt2var;
  }
  for (i = 0, s = sp.begin(); i < L; i++, s += nst) {
    checkUserInterrupt();
    ggmperm(p, nm, xx.begin());
    ggscale(p, n, m, true, score, iter, xx.begin(), zs.begin(), z.begin(),
            l.begin(), r.begin(), iw.begin(), w.begin());
    for (j = 0, zi = z.begin(); j < nm; j += n, zi += pn) {
      cstat(p, n, zi, st1.begin(), w.begin());
      for (k = 0; k < nst; k++)
        s[k] = std::max(s[k], st1[k]);
    }
  }
  ggscale(p, n, m, true, score, iter, x.begin(), zs.begin(), z.begin(),
          l.begin(), r.begin(), iw.begin(), w.begin());
  for (j = 0, zi = z.begin(), s = st.begin(); j < nm;
       j += n, zi += pn, s += nst) {
    cstat(p, n, zi, st1.begin(), w.begin());
    for (k = 0; k < nst; k++)
      s[k] = st1[k];
  }
  return List::create(_["center"] = l, _["scatter"] = r, _["statistic"] = st,
                      _["sp"] = sp);
}

// [[Rcpp::export]]
List ggscore2mshewhart(NumericVector x, std::string stat, int L = 1000) {
  IntegerVector dim = x.attr("dim");
  int i, j, k, p = dim[0], n = dim[1], m = dim(2), nm = n * m, pn = p * n,
               nst = (stat != "T2Var") ? 1 : 2;
  double *zi, *s;
  NumericVector xx = clone(x), st(nst * m), st1(2), sp(nst * L),
                w(6 * p + 2 * p * p);
  void (*cstat)(int, int, double *, double *, double *);
  if (stat != "T2Var") {
    cstat = (stat == "T2") ? &ggt2 : &ggvar;
  } else {
    st.attr("dim") = IntegerVector::create(2, m);
    sp.attr("dim") = IntegerVector::create(2, L);
    cstat = &ggt2var;
  }
  for (i = 0, s = sp.begin(); i < L; i++, s += nst) {
    checkUserInterrupt();
    ggmperm(p, nm, xx.begin());
    for (j = 0, zi = xx.begin(); j < nm; j += n, zi += pn) {
      cstat(p, n, zi, st1.begin(), w.begin());
      for (k = 0; k < nst; k++)
        s[k] = std::max(s[k], st1[k]);
    }
  }
  for (j = 0, zi = x.begin(), s = st.begin(); j < nm;
       j += n, zi += pn, s += nst) {
    cstat(p, n, zi, st1.begin(), w.begin());
    for (k = 0; k < nst; k++)
      s[k] = st1[k];
  }
  return List::create(_["statistic"] = st, _["sp"] = sp);
}

// [[Rcpp::export]]
List ggglrchart(NumericVector x, bool onlymean = false, int L = 1000) {
  IntegerVector dim = x.attr("dim");
  int i, p = dim[0], n = dim[1], m = dim(2), nm = n * m;
  double *gg;
  NumericVector xx = clone(x), g(m), w((m + 5) * (p + p * p));
  NumericMatrix gp(m, L);
  ggglrt(p, n, m, onlymean, x.begin(), g.begin(), w.begin());
  for (i = 0, gg = gp.begin(); i < L; i++, gg += m) {
    checkUserInterrupt();
    ggmperm(p, nm, xx.begin());
    ggglrt(p, n, m, onlymean, xx.begin(), gg, w.begin());
  }
  return List::create(_["glr"] = g, _["glr.perm"] = gp);
}
