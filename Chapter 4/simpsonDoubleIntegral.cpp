#include "temp.cpp"
using namespace std;

db simpson_double_integral(db a, db b, function<db(db)> c, function<db(db)> d, function<db(db, db)> f, int n, int m) {
    db h = (b - a) / db(n);
    db J1 = 0, J2 = 0, J3 = 0;
    for (int i = 0; i <= n; ++i) {
        db x = a + h * i;
        db HX = (d(x) - c(x)) / db(m);
        db K1 = f(x, c(x)) + f(x, d(x)), K2 = 0, K3 = 0;
        for (int j = 1; j < m; ++j) {
            db y = c(x) + HX * j;
            db Q = f(x, y);
            if (j % 2 == 0) K2 += Q;
            else K3 += Q;
        }
        db L = (K1 + K2 * 2 + K3 * 4) * HX / 3.0;
        if (i == 0 || i == n) J1 += L;
        else if (i % 2 == 0) J2 += L;
        else J3 += L;
    }
    return h * (J1 + 2 * J2 + 4 * J3) / 3.0;
}

int main() {
    db a = 0, b = 1;
    function<db(db)> c = [](db x) { return 0; };
    function<db(db)> d = [](db x) { return sqrt(1.0 - x * x); };
    function<db(db, db)> sigma = [](db x, db y) { return exp(-x * x - y * y); };
    function<db(db, db)> f1 = [&](db x, db y) { return x * sigma(x, y); };
    function<db(db, db)> f2 = [&](db x, db y) { return x * sigma(x, y); };
    int n = 14, m = 14;
    cout << fixed << setprecision(10) << simpson_double_integral(a, b, c, d, f1, n, m) / simpson_double_integral(a, b, c, d, sigma, n, m) << endl;
    cout << fixed << setprecision(10) << simpson_double_integral(a, b, c, d, f2, n, m) / simpson_double_integral(a, b, c, d, sigma, n, m) << endl;
    return 0;
}