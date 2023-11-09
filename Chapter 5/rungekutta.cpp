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
    db a = 1, b = 2;
    int N = 10;
    db initial_condition = -db(1)/log(2);
    auto f = [](db t, db y) {
        return y * y / (t + 1.0);
    };
    auto alphas = rungekutta(a, b, N, initial_condition, f);
    vd t(sz(alphas));
    F0R(i, sz(alphas)) t[i] = alphas[i].f;

    cout << endl << "Two Step Explicit Method" << endl;
    vd twoStep(11);
    F0R(i, 2) twoStep[i] = alphas[i].s;
    FOR(i, 1, 10) {
        twoStep[i + 1] = twoStep[i] + 0.1 / 2 * (3.0 * f(t[i], twoStep[i]) - f(t[i-1], twoStep[i-1]));
    }
    F0R(i, sz(twoStep)) cout << t[i] << ", " << twoStep[i] << endl;

    cout << endl << "Three Step Explicit Method" << endl;
    vd threeStep(11);
    F0R(i, 3) threeStep[i] = alphas[i].s;
    FOR(i, 2, 10) {
        threeStep[i + 1] = threeStep[i] + 0.1 / 12.0 * (23.0 * f(t[i], threeStep[i]) - 16.0 * f(t[i-1], threeStep[i-1]) + 5.0 * f(t[i-2], threeStep[i-2]));
    }
    F0R(i, sz(threeStep)) cout << t[i] << ", " << threeStep[i] << endl;

    cout << endl << "Four Step Explicit Method" << endl;
    vd fourStep(11);
    F0R(i, 4) fourStep[i] = alphas[i].s;
    FOR(i, 3, 10) {
        fourStep[i + 1] = fourStep[i] + 0.1 / 24.0 * (55.0 * f(t[i], fourStep[i]) - 59.0 * f(t[i-1], fourStep[i-1]) + 37.0 * f(t[i-2], fourStep[i-2]) - 9.0 * f(t[i-3], fourStep[i-3]));
    }
    F0R(i, sz(fourStep)) cout << t[i] << ", " << fourStep[i] << endl;

    cout << endl << "Five Step Explicit Method" << endl;
    vd fiveStep(11);
    F0R(i, 5) fiveStep[i] = alphas[i].s;
    FOR(i, 4, 10) {
        fiveStep[i + 1] = fiveStep[i] + 0.1 / 720.0 * (1901.0 * f(t[i], fiveStep[i]) - 2774.0 * f(t[i-1], fiveStep[i-1]) + 2616.0 * f(t[i-2], fiveStep[i-2]) - 1274.0 * f(t[i-3], fiveStep[i-3]) + 251.0 * f(t[i-4], fiveStep[i-4]));
    }
    F0R(i, sz(fiveStep)) cout << t[i] << ", " << fiveStep[i] << endl;

    cout << endl << "Real values" << endl;
    F0R(i, sz(t)) cout << t[i] << ", " << -1.0 / log(t[i] + 1) << endl;

    return 0;
}