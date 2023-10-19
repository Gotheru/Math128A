#include "temp.cpp"
using namespace std;

tcT> V<V<T>> Romberg_Integration(db a, db b, int n, function<T(db)> f) {
    V<V<T>> R(n, V<T>(n));
    db h = b - a;
    R[0][0] = h * (f(a) + f(b)) / 2.0;
    for (int i = 1; i < n; ++i) {
        R[i][0] = R[i - 1][0];
        for (int k = 1; k <= (1 << (i - 1)); ++k) R[i][0] += h * f(a + h * (db(k) - 0.5));
        R[i][0] /= 2.0;
        h /= 2.0;
        for (int j = 1; j <= i; ++j) R[i][j] = R[i][j - 1] + (R[i][j - 1] - R[i - 1][j - 1]) / db((1 << (2 * j)) - 1);
    }
    return R;
}

int main() {
    db a = 0, b = 1;
    function<db(db)> f = [](db x) { return 2.0 / sqrt(PI) * exp(-x * x); };
    int n = 5;
    auto R = Romberg_Integration(a, b, n, f);
    each(r, R) {
        each(ri,  r) cout << fixed << setprecision(10) << ri << ' ';
        cout << endl;
    }
    return 0;
}