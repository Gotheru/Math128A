#include "temp.cpp"
using namespace std;

// Algorithm 5.5
vector<pair<db, db>> adams(db a, db b, db initial_condition, db TOL, db hmax, db hmin, function<db(db,db)> f) {
    
    function<void(db,int,vd&,vd&)> RK4 = [&](db h, int i, vd& v, vd& x) {
        while (sz(v) < i + 4) v.eb();
        while (sz(x) < i + 4) x.eb();
        for (int j = i + 1; j <= i + 3; ++j) {
            db K1 = h * f(x[j-1], v[j-1]);
            db K2 = h * f(x[j-1] + h/2, v[j-1] + K1/2);
            db K3 = h * f(x[j-1] + h/2, v[j-1] + K2/2);
            db K4 = h * f(x[j-1] + h, v[j-1] + K3);
            v[j] = v[j-1] + (K1 + K2 * 2 + K3 * 2 + K4) / 6.0;
            x[j] = x[j-1] + h;
        }
    };
    V<pair<db,db>> ans;
    vd T(4), W(4);
    T[0] = a, W[0] = initial_condition;
    db h = hmax;
    int FLAG = 1, LAST = 0;
    ans.pb(mp(T[0],W[0]));
    RK4(h, 0, W, T);
    int NFLAG = 1;
    int i = 4;
    db t = T[3] + h;
    while (FLAG) {
        db WP = W[i-1] + h/24.0*(f(T[i-1],W[i-1])*55.0-f(T[i-2],W[i-2])*59.0+f(T[i-3],W[i-3])*37.0-f(T[i-4],W[i-4])*9);
        db WC = W[i-1] + h/24.0*(f(t,WP)*9.0+f(T[i-1],W[i-1])*19.0-f(T[i-2],W[i-2])*5.0+f(T[i-3],W[i-3]));
        db sigma = abs(WC-WP)*19.0/h/270.0;
        if (sigma <= TOL) { // Step 6
            // Step 7
            while (sz(W) <= i) W.eb();
            while (sz(T) <= i) T.eb();
            W[i] = WC;
            T[i] = t;
            // Step 8
            while (sz(ans) <= i) ans.eb();
            if (NFLAG) {
                for (int j = i - 3; j <= i; ++j) ans[j] = mp(T[j], W[j]);
            } else {
                ans[i] = mp(T[i], W[i]);
            }
            // Step 9
            if (LAST) break;
            // Step 10
            ++i;
            NFLAG = 0;
            // Step 11
            if (sigma <= 0.1 * TOL || T[i-1] + h > b) {
                // Step 12
                db q = pow(TOL / 2 / sigma, 0.25);
                // Step 13
                if (q > 4) h *= 4;
                else h *= q;
                // Step 14
                if (h > hmax) h = hmax;
                // Step 15
                if (T[i-1]+h*4 > b)
                    h = (b - T[i-1]) / 4, LAST = 1;
                // Step 16
                RK4(h, i - 1, W, T);
                NFLAG = 1;
                i = i + 3;
            }
        } else {
            // Step 17
            db q = pow(TOL / 2 / sigma, 0.25);
            // Step 18
            if (q < 0.1) h *= 0.1;
            else h *= q;
            // Step 19
            if (h < hmin) {

            } else if (NFLAG) {
                i = i-3;
                RK4(h, i - 1, W, T);
                i += 3;
                NFLAG = 1;
            }
        }
        t = T[i-1] + h;
    }
    return ans;
}

int main() {
    db a = 0, b = 0.5, initial_condition = 0;
    db TOL = 1E-4, hmin = 0.01, hmax = 0.1;
    function<db(db,db)> f = [](db t, db y) { return sin(t) + exp(-t); };
    each(p, adams(a, b, initial_condition, TOL, hmax, hmin, f)) {
        cout << p.f << ", " << p.s << endl;
    }
    return 0;
}