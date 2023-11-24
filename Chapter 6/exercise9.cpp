#include "temp.cpp"
#include "matrix.cpp"
#include "gauss.cpp"

Mat submat(Mat const& A, int i1, int i2, int j1, int j2) {
    Mat B = makeMat(i2 - i1, j2 - j1);
    FOR(i, i1, i2) FOR(j, j1, j2) B[i-i1][j-j1] = A[i][j];
    return B;
}

void printMatrix(Mat const& A) {
    each(a, A) cout << a << endl;
    cout << endl;
}

int main() {
    Mat A = {
        { 1, 2, -1 },
        { 3, -4, -3 },
        { 6, 5, 0 }
    };
    Mat B = { 
        { 2, -1, 7, 0 },
        { 3, 0, 4, 5 },
        { -2, 1, -3, 1 }
    };
    Mat A11 = submat(A, 0, 2, 0, 2);
    Mat A12 = submat(A, 0, 2, 2, 3);
    Mat A21 = submat(A, 2, 3, 0, 2);
    Mat A22 = submat(A, 2, 3, 2, 3);
    Mat B11 = submat(B, 0, 2, 0, 3);
    Mat B12 = submat(B, 0, 2, 3, 4);
    Mat B21 = submat(B, 2, 3, 0, 3);
    Mat B22 = submat(B, 2, 3, 3, 4);

    printMatrix(A * B);
    printMatrix(A11 * B11 + A12 * B21);
    printMatrix(A11 * B12 + A12 * B22);
    printMatrix(A21 * B11 + A22 * B21);
    printMatrix(A21 * B12 + A22 * B22);

    return 0;
}