#include "temp.cpp"
using namespace std;

vector<pair<db, db>> rungekuttafehlberg(db a, db b, db initial_condition, db TOL, db hmax, db hmin, function<db(db, db)> f) {
    vector<pair<db, db>> ans;
    db t = a, w = initial_condition, h = hmax;
    int FLAG = 1;
    ans.pb(mp(t,w));
    while (FLAG) {
        db K1 = h * f(t, w);
        db K2 = h * f(t + h / 4, w + K1 / 4);
        db K3 = h * f(t + h * 3 / 8, w + K1 * 3 / 32 + K2 * 9 / 32);
        db K4 = h * f(t + h * 12 / 13, w + K1 * 1932 / 2197 - K2 * 7200 / 2197 + K3 * 7296 / 2197);
        db K5 = h * f(t + h, w + K1 * 439 / 216 - K2 * 8 + K3 * 3680 / 513 - K4 * 845 / 4104);
        db K6 = h * f(t + h / 2, w - K1 *8 / 27 + K2 * 2 - K3 * 3544 / 2565 + K4 * 1859 / 4104 - K5 * 11 / 40);
        db R = 1.0 / h * abs(K1 / 360 - K3 * 128 / 4275 - K4 * 2197 / 75240 + K5 / 50 + K6 * 2 / 55);
        if (R <= TOL) {
            t = t + h;
            w = w + K1 * 25 / 216 + K3 * 1408 / 2565 + K4 * 2197 / 4104 - K5 / 5;
            ans.pb(mp(t,w));
        }
        else {
            db delta = 0.84 * pow(TOL / R, 0.25);
            h = min(min(db(4), max(delta, db(0.1))) * h, hmax);
        }
        if (t >= b) FLAG = 0;
        else if (t + h > b) h = b - t;
        // else if (h < hmin) return {};
    }
    return ans;
}

int main() {
    db a = 1, b = 4, initial_condition = -1.0/log(2);
    db TOL = 1E-6, hmin = 0.05, hmax = 0.5;
    function<db(db,db)> f = [](db t, db y) { return y * y / (t + 1); };
    each(p, rungekuttafehlberg(a, b, initial_condition, TOL, hmin, hmax, f)) {
        cout << p.f << ", " << p.s << endl;
    }
    return 0;
}