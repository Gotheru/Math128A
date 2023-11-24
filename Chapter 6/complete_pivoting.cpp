#include "temp.cpp"
using namespace std;

template <class T>
vector<T> partial_pivoting(int n, vector<vector<T>> a) {
    assert(sz(a) == n);
    assert(sz(a[0]) == n + 1);
    vi nrow(n), ncolumn(n);
    iota(all(nrow), 0);
    iota(all(ncolumn), 0);
    cout << "nrow = " << nrow << endl;
    cout << "ncolumn = " << ncolumn << endl;
    for (int k = 0; k < n - 1; ++k) {
        int pivot_i = k, pivot_j = k;
        for (int i = k; i < n; ++i) {
            for (int j = k; j < n; ++j) {
                if (abs(a[nrow[pivot_i]][ncolumn[pivot_j]]) < abs(a[nrow[i]][ncolumn[j]])) {
                    pivot_i = i, pivot_j = j;
                }
            }
        }
        swap(nrow[k], nrow[pivot_i]);
        swap(ncolumn[k], ncolumn[pivot_j]);
        for (int j = k + 1; j < n; ++j) {
            auto m = a[nrow[j]][ncolumn[k]] / a[nrow[k]][ncolumn[k]];
            a[nrow[j]] = a[nrow[j]] - m * a[nrow[k]];
        }
        cout << "nrow = " << nrow << endl;
        cout << "ncolumn = " << ncolumn << endl;
    }
    if (a[nrow[n-1]][n-1] == 0) {
        cout << "No unique solution exists" << endl;
        return {};
    }
    vector<T> x(n);
    x[n-1] = a[nrow[n-1]][n] / a[nrow[n-1]][ncolumn[n-1]];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = a[nrow[i]][n];
        for (int j = i + 1; j < n; ++j) x[i] -= a[nrow[i]][ncolumn[j]] * x[j];
        x[i] /= a[nrow[i]][ncolumn[i]];
    }
    vector<T> ans(n);
    for (int i = 0; i < n; ++i) ans[ncolumn[i]] = x[i];
    return ans;
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