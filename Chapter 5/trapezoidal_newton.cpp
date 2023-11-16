#include "temp.cpp"
using namespace std;

V<pair<db, db>> trapezoidal_newton(db a, db b, int N, db initial_condition, db TOL, int M, function<db(db, db)> f, function<db(db, db)> fy) {
    vector<pair<db, db>> ans;
    db h = (b - a) / db(N);
    db t = a;
    db w = initial_condition;
    ans.pb(mp(t, w));
    F0R(i, N) {
        db k1 = w + h / 2 * f(t, w);
        db w0 = k1;
        int j = 1, FLAG = 0;
        while (!FLAG) {
            w = w0 - (w0 - h / 2 * f(t + h, w0) - k1) / (db(1) - h / 2 * fy(t + h, w0));
            if (abs(w - w0) < TOL) FLAG = 1;
            else {
                ++j;
                if (j > M) {
                    cout << "maximum number of iterations exceeded" << endl;
                }
            }
        }
        t += h;
        ans.pb(mp(t, w));
    }
    return ans;
}

int main() {
    db a = 0, b = 1, initial_condition = exp(1);
    int N = 10;
    db TOL = 1E-3;
    int M = 100;
    function<db(db, db)> f = [](db t, db y) { return -y * 10 + t * 10 + 1; };
    function<db(db, db)> fy = [](db t, db y) { return -10; };
    auto ans = trapezoidal_newton(a, b, N, initial_condition, TOL, M, f, fy);
    function<db(db)> y = [](db t) { return exp(-t * 10 + 1) + t; };
    each(row, ans) {
        cout << row.f << ' ' << row.s << ' ' << y(row.f) << '\n';
    }
    return 0;
}