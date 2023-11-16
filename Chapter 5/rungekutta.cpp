#include "temp.cpp"
using namespace std;

// Algorithm 5.2
vector<pair<db, db>> rungekutta(db a, db b, int N, db initial_condition, function<db(db, db)> f) {
    vector<pair<db, db>> ans;
    db h = (b - a) / db(N), t = a, w = initial_condition;
    ans.pb(mp(t,w));
    F0R(i, N) {
        db K1 = h * f(t, w);
        db K2 = h * f(t + h / 2, w + K1 / 2);
        db K3 = h * f(t + h / 2, w + K2 / 2);
        db K4 = h * f(t + h, w + K3);
        w = w + (K1 + K2 * 2.0 + K3 * 2.0 + K4) / 6.0;
        t += h;
        ans.pb(mp(t,w));
    }
    return ans;
}

int main() {
    db a = 0, b = 1; int N = 10;
    db initial_condition = exp(1);
    function<db(db, db)> f = [](db t, db y) { return -y * 10 + t * 10 + db(1); };
    function<db(db)> y = [](db t) { return exp(-10 * t + 1) + t; };
    auto ans = rungekutta(a, b, N, initial_condition, f);
    each(row, ans) {
        cout << fixed << setprecision(5) << row.f << ' ' << row.s << ' ' << y(row.f) << endl;
    }
    return 0;
}