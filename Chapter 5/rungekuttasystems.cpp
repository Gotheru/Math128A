#include "temp.cpp"
using namespace std;

vd operator + (vd v, vd w) {
    assert(sz(v) == sz(w));
    F0R(i, sz(v)) v[i] += w[i];
    return v;
}

vd operator * (db c, vd v) {
    F0R(i, sz(v)) v[i] *= c;
    return v;
}

V<pair<db, vd>> runge_kutta_systems(db a, db b, int m, int N, vd const& initial_conditions, V<function<db(db, vd)>> const& f) {
    db h = (b - a) / db(N);
    db t = a;
    V<pair<db, vd>> ans;
    vd w(m);
    F0R(j, m) w[j] = initial_conditions[j];
    ans.pb(mp(t, w));
    for (int i = 1; i <= N; ++i) {
        V<vd> k(4, vd(m));
        F0R(j, m) k[0][j] = h * f[j](t, w);
        F0R(j, m) k[1][j] = h * f[j](t + h / 2, w + db(0.5) * k[0]);
        F0R(j, m) k[2][j] = h * f[j](t + h / 2, w + db(0.5) * k[1]);
        F0R(j, m) k[3][j] = h * f[j](t + h, w + k[2]);
        F0R(j, m) w[j] += (k[0][j] + k[1][j] * 2 + k[2][j] * 2 + k[3][j]) / db(6);
        t += h;
        ans.pb(mp(t, w));
    }
    return ans;
}

ostream& operator << (ostream& o, vd const& v) {
    o << "( ";
    F0R(i, sz(v) - 1) o << v[i] << ", ";
    if (sz(v)) o << v.back() << ' ';
    return o << ')';
}

int main() {
    // 2b
    /*
    db a = 0, b = 2;
    int m = 2;
    int N = 10;
    vd initial_conditions = {-3, 5};
    function<db(db, vd)> f0 = [](db t, vd u) {
        return u[0] / 9 - u[1] * 2 / 3 - t * t / 9 + db(2) / 3;
    };
    function<db(db, vd)> f1 = [](db t, vd u) {
        return u[1] + t * 3 - db(4);
    };
    V<function<db(db, vd)>> f = { f0, f1 };
    auto ans = runge_kutta_systems(a, b, m, N, initial_conditions, f);
    function<db(db)> u1 = [](db t) { return - db(3) * exp(t) + t * t; };
    function<db(db)> u2 = [](db t) { return db(4) * exp(t) - 3 * t + 1; };
    function<vd(db)> u = [&](db t) { return vd({u1(t), u2(t)}); };
    each(row, ans) {
        cout << fixed << setprecision(5) << row.f << ' ' << row.s << ' ' << u(row.f) << endl;
    }
    */

    // 4B
    /*
    db a = 1, b = 3;
    int m = 2;
    int N = 10;
    vd initial_conditions = {4, 3};
    function<db(db, vd)> f0 = [](db t, vd u) {
        return u[1];
    };
    function<db(db, vd)> f1 = [](db t, vd u) {
        return (-t * 3 + u[0] * 4 - t * u[1]) / t / t;
    };
    V<function<db(db, vd)>> f = { f0, f1 };
    auto ans = runge_kutta_systems(a, b, m, N, initial_conditions, f);
    function<db(db)> y = [](db t) { return t * t * 2 + t + db(1) / t / t; };
    each(row, ans) {
        cout << fixed << setprecision(5) << row.f << ' ' << row.s[0] << ' ' << y(row.f) << endl;
    }
    */
    // 5
    /*
    db a = 0, b = 4;
    int m = 2;
    int N = 40;
    vd initial_conditions = {1000, 500};
    db k1 = 3, k2 = 0.002, k3 = 0.0006, k4 = 0.5;
    function<db(db, vd)> f0 = [&](db t, vd x) {
        return k1 * x[0] - k2 * x[0] * x[1];
    };
    function<db(db, vd)> f1 = [&](db t, vd x) {
        return k3 * x[0] * x[1] - k4 * x[1];
    };
    V<function<db(db, vd)>> f = { f0, f1 };
    auto ans = runge_kutta_systems(a, b, m, N, initial_conditions, f);
    each(row, ans) {
        cout << fixed << setprecision(5) << row.f << ", " << row.s[0] << endl;
    }

    cout << endl << endl;

    each(row, ans) {
        cout << fixed << setprecision(5) << row.f << ", " << row.s[1] << endl;
    }
    */
    // 6
    db a = 0, b = 4;
    int m = 2;
    int N = 100;
    vd initial_conditions = {10000, 10000};
    function<db(db, vd)> f0 = [&](db t, vd x) {
        return x[0] * (db(4) - x[0] * 0.0003 - x[1] * 0.0004);
    };
    function<db(db, vd)> f1 = [&](db t, vd x) {
        return x[1] * (db(2) - x[0] * 0.0002 - x[1] * 0.0001);
    };
    V<function<db(db, vd)>> f = { f0, f1 };
    auto ans = runge_kutta_systems(a, b, m, N, initial_conditions, f);
    each(row, ans) {
        cout << fixed << setprecision(5) << row.f << ", " << row.s[0] << endl;
    }

    cout << endl << endl;

    each(row, ans) {
        cout << fixed << setprecision(5) << row.f << ", " << row.s[1] << endl;
    }

    return 0;
}