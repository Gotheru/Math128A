#include "temp.cpp"
using namespace std;

template <class T>
vector<T> partial_pivoting(int n, vector<vector<T>> a) {
    assert(sz(a) == n);
    assert(sz(a[0]) == n + 1);
    vi nrow(n);
    iota(all(nrow), 0);
    vector<T> s(n);
    for (int i = 0; i < n; ++i) {
        s[i] = 0;
        for (int j = 0; j < n; ++j) {
            s[i] = max(a[i][j], s[i]);
        }
        if (s[i] == 0) {
            cout << "No unique solution exists\n";
            return {};
        }
    }
    cout << "nrow = " << nrow << endl;
    for (int i = 0; i < n - 1; ++i) {
        int p = i;
        for (int j = i; j < n; ++j) {
            if (abs(a[nrow[p]][i]) / s[nrow[p]] < abs(a[nrow[j]][i]) / s[nrow[p]]) p = j;
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
    int n = 3;
    V<V<db>> a = {
        { 1, 1, -1, 0 },
        { 0, 12, -1, 4 },
        { 2, 1, 1, 5 }
    };
    cout << partial_pivoting(n, a) << endl;
}