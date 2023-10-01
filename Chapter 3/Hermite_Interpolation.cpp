#include "temp.cpp"
#include "Polynomial.cpp"

tcT> Polynomial<T> hermite_interpolation(V<T> const& x, V<T> const& f, V<T> const& df) {

    int N = sz(x);

    assert(N != 0);
    assert(N == sz(f));
    assert(N == sz(df));

    V<V<T>> Q(2 * N, V<T>(2 * N, 0));
    V<T> z(2 * N);
    F0R(i, N) {
        z[2 * i] = x[i];
        z[2 * i + 1] = x[i];
        Q[2 * i][0] = Q[2 * i + 1][0] = f[i];
        Q[2 * i + 1][1] = df[i];
        if (i != 0) {
            Q[2 * i][1] = (Q[2 * i][0] - Q[2 * i - 1][0]) / (z[2 * i] - z[2 * i - 1]);
        }
    }

    FOR(i, 2, 2 * N) {
        FOR(j, 2, i + 1) {
            Q[i][j] = (Q[i][j - 1] - Q[i - 1][j - 1]) / (z[i] - z[i - j]);
        }
    }

    V<T> ans;
    F0R(i, 2 * N) ans.pb(Q[i][i]);

    return ans;

}

db f(db x) { return exp(0.1 * x * x); }
db df(db x) { return 0.2 * x * exp(0.1 * x * x); }

int main() {
	setIO();
    V<db> x = { 1, 1.5 };
    V<db> fx; each(xi, x) fx.pb(f(xi));
    V<db> dfx; each(xi, x) dfx.pb(df(xi));
    Polynomial H = hermite_interpolation(x, fx, dfx);
    each(h, H.coeffs) cout << setprecision(4) << h << ' ';
    cout << endl;
    V<db> z(2 * sz(x));
    F0R(i, sz(x)) z[2 * i] = z[2 * i + 1] = x[i];
    cout << setprecision(10) << H.evaluate_difference(z, 1.25) << '\n';
    cout << setprecision(10) << f(1.25) - H.evaluate_difference(z, 1.25) << '\n';
    return 0;
}