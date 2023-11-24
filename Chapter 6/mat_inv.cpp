#include "gauss.cpp"

int main() {
    Mat A = {
        { 2, 0, 1, 2 },
        { 1, 1, 0, 2 },
        { 2, -1, 3, 1 },
        { 3, -1, 4, 3 }
    };
    cout << "Determinant = " << gauss(A).f << endl;
    Mat A_1 = inv(A);
    each(a, A) cout << a << endl;
    return 0;
}