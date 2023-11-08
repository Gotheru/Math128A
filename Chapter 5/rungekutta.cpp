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
        db w = w + (K1 + K2 * 2 + K3 * 3 + K4) / 6.0;
        t += h;
        ans.pb(mp(t,w));
    }
    return ans;
}