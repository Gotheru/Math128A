#include "temp.cpp"

using T = db;
vector<T> substitution(int n, vector<vector<T>> a) {
    assert(sz(a) == n);
    assert(sz(a[0]) == n + 1);
    for (int i = 0; i < n - 1; ++i) {
        each(ai, a) cout << ai << endl;
        cout << endl;
        int p = i;
        while (p < n && a[p][i] == 0) ++p;
        cout << p << endl;
        if (p == n) {
            cout << "No unique solution exists\n";
            return {};
        }
        if (i != p) {
            cout << "Swapped rows i  = " << i << " and p = " << p << '\n';
            a[i].swap(a[p]);
        }
        for (int j = i + 1; j < n; ++j) {
            T m = a[j][i] / a[i][i];
            a[j] = a[j] - m * a[i];
        }
    }
    each(ai, a) cout << ai << endl;
    cout << endl;
    if (a[n - 1][n - 1] == 0) {
        cout << "No unique solution exists\n";
        return {};
    }
    vector<T> x(n);
    x[n-1] = a[n-1][n] / a[n-1][n-1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = a[i][n];
        for (int j = i + 1; j < n; ++j) x[i] -= x[j] * a[i][j];
        x[i] /= a[i][i];
    }
    return x;
}

int main() {
    int n = 5;
    db h = 0.25;
    vector<vd> a(5, vd(6));
    vd x = { 0, 0.25, 0.5, 0.75, 1 };
    auto f = [](db x) { return x * x; };
    auto K = [](db x, db t) { return exp(abs(x- t)); };
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i][j] = K(x[i], x[j]);
            if (j != 0 || j != n - 1) a[i][j] *= 2;
            a[i][j] *= h / 2;
        }
        a[i][i] -= 1;
        a[i][n] = -f(x[i]);
    }
    cout << substitution(n, a) << endl;
}