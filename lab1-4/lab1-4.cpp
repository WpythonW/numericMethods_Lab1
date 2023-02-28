#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

int n;
double eps;

struct pair_answ {
    int i, j;
    double max, sqr;
};

void print_matr(double** M) {
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) {
            cout << M[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

double** transpose(double** M) {
    double** M_out = new double* [n];

    for (int i = 0; i < n; i++)
    {
        M_out[i] = new double[n];
        for (int j = 0; j < n; j++) {
            M_out[i][j] = M[j][i];
        }
    }
    return M_out;
}

double** multiple(double** M1, double** M2) {
    double** res = new double* [n];

    for (int i = 0; i < n; ++i) {
        res[i] = new double[n];
        for (int j = 0; j < n; ++j) {
            res[i][j] = 0;
            for (int k = 0; k < n; ++k)
                res[i][j] += M1[i][k] * M2[k][j];
        }
    }
    return res;
}

struct pair_answ find_max(double **A) {
    double max = 0, summ = 0;
    pair_answ answ;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (abs(A[i][j]) > abs(max)) {
                max = A[i][j];
                answ.i = i;
                answ.j = j;
                summ += A[i][j] * A[i][j];
            }
        }
    }
    answ.max = max;
    answ.sqr = sqrt(summ);
    return answ;
}

double** findU(int mi, int mj, double phi) {
    double** U = new double* [n];
    for (int i = 0; i < n; i++) {
        U[i] = new double[n];
        for (int j = 0; j < n; j++)
            if (i != j) U[i][j] = 0;
        U[i][i] = 1;
        if (i == mi) {
            U[mi][mi] = cos(phi);
            U[mi][mj] = -sin(phi);
        }
        else if (i == mj){
            U[mj][mj] = cos(phi);
            U[mj][mi] = sin(phi);
        }
    }
    return U;
}

void rotation_method(double** A) {
    double phi;
    double** U;
    double** X = new double*[n];


    pair_answ answ = find_max(A);
    if (A[answ.i][answ.i] - A[answ.j][answ.j] != 0)
        phi = 0.5 * atan(2 * answ.max / (A[answ.i][answ.i] - A[answ.j][answ.j]));
    else
        phi = M_PI / 4;
    U = findU(answ.i, answ.j, phi);

    copy(U, U + n, X);

    A = multiple(multiple(transpose(U), A), U);

    answ = find_max(A);



    while (answ.sqr > eps) {
        if (A[answ.i][answ.i] - A[answ.j][answ.j] != 0)
            phi = 0.5 * atan(2 * answ.max / (A[answ.i][answ.i] - A[answ.j][answ.j]));
        else
            phi = M_PI / 4;

        U = findU(answ.i, answ.j, phi);
        X = multiple(X, U);
        A = multiple(multiple(transpose(U), A), U);

        answ = find_max(A);
    } 

    print_matr(A);
    cout << "\n";
    print_matr(X);

}

int main() {
    ifstream inFile;

    inFile.open("matrix.txt");
    inFile >> eps;
    inFile >> n;
    double** A = new double* [n];
    for (int i = 0; i < n; i++)
        A[i] = new double[n];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
        {
            inFile >> A[i][j];
            cout << A[i][j] << " ";
        }
        cout << "\n";
    }

    rotation_method(A);
}
