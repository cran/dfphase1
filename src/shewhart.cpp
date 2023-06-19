#define STRICT_R_HEADERS
#include <Rcpp.h>
using namespace Rcpp;

namespace {

void ggperm(int n, double *y) {
    int i;
    double yi;
    while (n > 1) {
        i = floor(n * unif_rand());
        n -= 1;
        yi = y[i];
        y[i] = y[n];
        y[n] = yi;
    }
}

inline double ggmedian(int n, double *x, double *w) {
    int half = n / 2;
    std::copy(x, x + n, w);
    std::nth_element(w, w + half, w + n);
    if (n == 2 * half) {
        return (*std::max_element(w, w + half) + w[half]) / 2;
    } else {
        return w[half];
    }
}

struct Comparator {
    Comparator(const double *x) : data(x) {}
    bool operator()(int left, int right) const {
        return data[left] < data[right];
    }
    const double *data;
};

inline void rank(int n, double *x, double *r, int *a) {
    int i, j, s;
    double d;
    for (i = 0; i < n; i++) a[i] = i;
    std::sort(a, a + n, Comparator(x));
    for (i = 0; i < n;) {
        for (s = i, j = i + 1; (j < n) && (x[a[j]] <= x[a[i]]); j++) s += j;
        for (d = 1 + ((double)s) / (j - i); i < j; i++) r[a[i]] = d;
    }
}

inline void ggcolmeans(NumericMatrix x, NumericVector xbar) {
    int n = x.nrow(), m = x.ncol(), i, j;
    double a, *xi = x.begin();
    for (i = 0; i < m; i++, xi += n) {
        for (j = 0, a = 0.0; j < n; j++) a += xi[j];
        xbar[i] = a / n;
    }
}

inline void ggcolxbars(NumericMatrix x, NumericVector xbar, NumericVector s) {
    int n = x.nrow(), m = x.ncol(), i, j;
    double a, b, c,
        *xi = x.begin(),
        c4 = M_SQRT2 * exp(R::lgammafn(0.5 * n) - R::lgammafn(0.5 * (n - 1)));
    for (i = 0; i < m; i++, xi += n) {
        for (j = 0, a = b = 0.0; j < n; j++) {
            a += c = xi[j];
            b += c * c;
        }
        xbar[i] = a / n;
        s[i] = sqrt(b - a * a / n) / c4;
    }
}

// length(w)=m if !aggr_with_mean
inline void horsexbars(NumericMatrix x, bool aggr_with_mean, NumericVector xbar,
                       NumericVector s, NumericVector est, NumericVector w) {
    ggcolxbars(x, xbar, s);
    if (aggr_with_mean) {
        est[0] = mean(xbar);
        est[1] = mean(s);
    } else {
        est[0] = ggmedian(xbar.size(), xbar.begin(), w.begin());
        est[1] = ggmedian(s.size(), s.begin(), w.begin());
    }
}

inline void horserank(NumericMatrix x, NumericVector l, NumericVector s,
                      NumericMatrix r, IntegerVector a) {
    int i, j, n = x.nrow(), m = x.ncol(), N = m * n, MN = m * (n - 1);
    double aq, bq, cq, mq, sq,
        med = ggmedian(N, x.begin(), r.begin()),
        c4 = exp(R::lgammafn(0.5 * MN) - R::lgammafn(0.5 * (MN - 1))) *
             sqrt(2.0 / (MN - 1.0));
    rank(N, x.begin(), r.begin(), a.begin());
    ggcolmeans(r, l);
    l = (l - (N + 1) / 2.0) / sqrt((m - 1) * (N + 1) / 12.0);
    for (i = 0; i < N; i++) r[i] = fabs(x[i] - med);
    rank(N, r.begin(), r.begin(), a.begin());
    for (i = 0, mq = sq = 0.0; i < m; i++) {
        for (j = 0, aq = bq = 0; j < n; j++) {
            aq += cq = r(j, i) * r(j, i);
            bq += cq * cq;
        }
        aq /= n;
        s[i] = aq;
        mq += aq;
        sq += (bq - n * aq * aq) / (n - 1);
    }
    mq /= m;
    sq = sqrt(sq / m) / c4;
    s = sqrt(static_cast<double>(n)) * (s - mq) / sq;
}

inline void horselepage(NumericMatrix x, NumericVector l, NumericVector s,
                        NumericVector lp, NumericMatrix r, IntegerVector a) {
    int i, j, n = x.nrow(), m = x.ncol(), N = m * n;
    double li, si,
        mid = 0.5 * (N + 1.0), w2 = N * (N + 1.0) / 12.0,
        mu = (N % 2) ? (N + 2.0) / 4.0 : (N + 1.0) * (N + 1.0) / (4.0 * N),
        v2 = (N % 2) ? N * (N * N - 4.0) / (48.0 * (N - 1.0))
                     : (N + 1.0) * (N * N + 3.0) / (48.0 * N);
    rank(N, x.begin(), r.begin(), a.begin());
    for (i = 0; i < m; i++) {
        for (j = 0, li = 0, si = 0; j < n; j++) {
            li += r(j, i);
            si += mid - fabs(r(j, i) - mid);
        }
        li /= n;
        si /= n;
        l[i] = n * (li - mid) * (li - mid) / w2;
        s[i] = n * (si - mu) * (si - mu) / v2;
        lp[i] = l[i] + s[i];
    }
}

inline void horsecucconi(NumericMatrix x, NumericVector l, NumericVector s,
                         NumericVector lp, NumericMatrix r, IntegerVector a) {
    int i, j, n = x.nrow(), m = x.ncol(), N = m * n;
    double li, si,
        mu = n * (N + 1.0) * (2.0 * N + 1 - 0) / 6.0,
        sigma =
            sqrt(static_cast<double>(n * (N - n) * (N + 1.0) * (2.0 * N + 1.0) *
                                     (8.0 * N + 11.0) / 180.0)),
        rho = -(30.0 * N + 14.0 * N * N + 19) /
              ((8.0 * N + 11) * (2.0 * N + 1.0)),
        N1 = N + 1.0, d1 = 2.0 * m, d2 = d1 * (1 - rho * rho);
    rank(N, x.begin(), r.begin(), a.begin());
    for (i = 0; i < m; i++) {
        for (j = 0, li = 0, si = 0; j < n; j++) {
            li += r(j, i) * r(j, i);
            si += (N1 - r(j, i)) * (N1 - r(j, i));
        }
        li = (li - mu) / sigma;
        si = (si - mu) / sigma - rho * li;
        l[i] = li * li / d1;
        s[i] = si * si / d2;
        lp[i] = l[i] + s[i];
    }
}

}  // namespace

// [[Rcpp::export]]
List ggxbars(NumericMatrix x, bool aggr_with_mean, int L) {
    int i, n = x.nrow(), m = x.ncol(), nm = n * m;
    double sn = sqrt(static_cast<double>(n));
    NumericVector xb(m), s(m), est(2), w(aggr_with_mean ? 0 : m);
    NumericMatrix xx = clone(x), pstat(3, L);
    for (i = 0; i < L; i++) {
        checkUserInterrupt();
        ggperm(nm, xx.begin());
        horsexbars(xx, aggr_with_mean, xb, s, est, w);
        xb = abs(xb - est[0]);
        pstat(0, i) = sn * (*std::max_element(xb.begin(), xb.end())) / est[1];
        pstat(1, i) = -(*std::min_element(s.begin(), s.end())) / est[1];
        pstat(2, i) = (*std::max_element(s.begin(), s.end())) / est[1];
    }
    horsexbars(x, aggr_with_mean, xb, s, est, w);
    return List::create(_["Xbar"] = xb, _["S"] = s, _["center"] = est[0],
                        _["scale"] = est[1], _["perm"] = pstat);
}

// [[Rcpp::export]]
NumericMatrix ggxbarsall(int n, int m, bool aggr_with_mean, int rep) {
    int i;
    double sn = sqrt(static_cast<double>(n));
    NumericVector xb(m), s(m), est(2), w(aggr_with_mean ? 0 : m);
    NumericMatrix x(n, m), stat(3, rep);
    for (i = 0; i < rep; i++) {
        checkUserInterrupt();
        std::generate(x.begin(), x.end(), norm_rand);
        horsexbars(x, aggr_with_mean, xb, s, est, w);
        xb = abs(xb - est[0]);
        stat(0, i) = sn * (*std::max_element(xb.begin(), xb.end())) / est[1];
        stat(1, i) = -(*std::min_element(s.begin(), s.end())) / est[1];
        stat(2, i) = (*std::max_element(s.begin(), s.end())) / est[1];
    }
    return stat;
}

// [[Rcpp::export]]
List ggrank(NumericMatrix x, int L) {
    int i, j, n = x.nrow(), m = x.ncol(), N = m * n;
    double A, B;
    NumericMatrix xx = clone(x), r(n, m), stat(2, L);
    IntegerVector a(N);
    NumericVector l(m), s(m);
    for (i = 0; i < L; i++) {
        checkUserInterrupt();
        ggperm(N, xx.begin());
        horserank(xx, l, s, r, a);
        for (j = 0, A = 0, B = 0; j < m; j++) {
            A = std::max(A, fabs(l[j]));
            B = std::max(B, fabs(s[j]));
        }
        stat(0, i) = A;
        stat(1, i) = B;
    }
    horserank(x, l, s, r, a);
    return List::create(_["lRank"] = l, _["sRank"] = s, _["perm"] = stat);
}

// [[Rcpp::export]]
NumericMatrix ggrankall(int n, int m, int rep) {
    int i, j;
    double A, B;
    NumericVector l(m), s(m);
    NumericMatrix x(n, m), r(n, m), stat(2, rep);
    IntegerVector a(m * n);
    for (i = 0; i < rep; i++) {
        checkUserInterrupt();
        std::generate(x.begin(), x.end(), unif_rand);
        horserank(x, l, s, r, a);
        for (j = 0, A = 0, B = 0; j < m; j++) {
            A = std::max(A, fabs(l[j]));
            B = std::max(B, fabs(s[j]));
        }
        stat(0, i) = A;
        stat(1, i) = B;
    }
    return stat;
}

// [[Rcpp::export]]
List gglepagecucconi(NumericMatrix x, int L, bool lepage) {
    int i, n = x.nrow(), m = x.ncol(), N = m * n;
    NumericMatrix xx = clone(x), r(n, m);
    IntegerVector a(N);
    NumericVector l(m), s(m), lp(m), stat(L);
    for (i = 0; i < L; i++) {
        checkUserInterrupt();
        ggperm(N, xx.begin());
        if (lepage)
            horselepage(xx, l, s, lp, r, a);
        else
            horsecucconi(xx, l, s, lp, r, a);
        stat[i] = *std::max_element(lp.begin(), lp.end());
    }
    if (lepage)
        horselepage(x, l, s, lp, r, a);
    else
        horsecucconi(x, l, s, lp, r, a);
    return List::create(_["Location"] = l, _["Scale"] = s, _["CS"] = lp,
                        _["perm"] = stat);
}

// [[Rcpp::export]]
NumericVector gglepagecucconiall(int n, int m, int rep, bool lepage) {
    int i;
    NumericVector l(m), s(m), lp(m), stat(rep);
    NumericMatrix x(n, m), r(n, m);
    IntegerVector a(m * n);
    for (i = 0; i < rep; i++) {
        checkUserInterrupt();
        std::generate(x.begin(), x.end(), unif_rand);
        if (lepage)
            horselepage(x, l, s, lp, r, a);
        else
            horsecucconi(x, l, s, lp, r, a);
        stat[i] = *std::max_element(lp.begin(), lp.end());
    }
    return stat;
}
