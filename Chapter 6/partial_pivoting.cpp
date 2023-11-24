#include "temp.cpp"
using namespace std;

template <class T>
vector<T> partial_pivoting(int n, vector<vector<T>> a) {
    assert(sz(a) == n);
    assert(sz(a[0]) == n + 1);
    vi nrow(n);
    iota(all(nrow), 0);
    cout << "nrow = " << nrow << endl;
    for (int i = 0; i < n - 1; ++i) {
        int p = i;
        for (int j = i; j < n; ++j) {
            if (abs(a[nrow[p]][i]) < abs(a[nrow[j]][i])) p = j;
        }
        if (nrow[i] != nrow[p]) {
            swap(nrow[i], nrow[p]);
        }
        for (int j = i + 1; j < n; ++j) {
            auto m = a[nrow[j]][i] / a[nrow[i]][i];
            a[nrow[j]] = a[nrow[j]] - m * a[nrow[i]];
        }
        cout << "nrow = " << nrow << endl;
    }
    if (a[nrow[n-1]][n-1] == 0) {
        cout << "No unique solution exists" << endl;
        return {};
    }
    vector<T> x(n);
    x[n-1] = a[nrow[n-1]][n] / a[nrow[n-1]][n-1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = a[nrow[i]][n];
        for (int j = i + 1; j < n; ++j) x[i] -= a[nrow[i]][j] * x[j];
        x[i] /= a[nrow[i]][i];
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
            if (j != 0 || j != n - 1) {
                if (j % 2 == 1) a[i][j] *= 4;
                else a[i][j] *= 2;
            }
            a[i][j] *= h / 3;
        }
        a[i][i] -= 1;
        a[i][n] = -f(x[i]);
    }
    cout << partial_pivoting(n, a) << endl;
}