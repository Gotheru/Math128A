#include "temp.cpp"
using namespace std;

const db EPS = 1e-10;

tcT> T S(db a, db b, function<T(db)> f) {
    return (b - a) / 6.0 * (f(a) + 4 * f((b + a) / 2.0) + f(b));
}

tcT> T adaptive_quadrature(db a, db b, function<T(db)> f, db TOL) {
    if (a + EPS > b) return 0;
    db m = (a + b) / 2.0;
    if (abs(S(a, b, f) - S(a, m, f) - S(m, b, f)) < 15 * TOL) return S(a, m, f) + S(m, b, f);
    else return adaptive_quadrature(a, m, f, TOL / 2) + adaptive_quadrature(m, b, f, TOL / 2);
}

int main() {
    db a = 0, b = 2 * PI;
    function<db(db)> f = [](db x) {
        return (cos(2 * x) - cos(3 * x)) / 5.0;
    };
    db TOL = 1e-5;
    cout << setprecision(12) << adaptive_quadrature(a, b, f, TOL) << endl;
    return 0;
}