#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

int n;
double eps;

struct complex
{
    double real;
    double im;
};

double sgn(double x)
{
    double sgn;
    if (x > 0.0L) sgn = 1.0L;
    if (x < 0.0L) sgn = -1.0L;
    if (x == 0.0L) sgn = 0.0L;
    return sgn;
}

void print(double** M) {
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) {
            cout << M[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
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

void ones(double** M) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            M[i][j] = 0;
        }
        M[i][i] = 1;
    }
}

double find_norm(double** M, int elem) {
    double summ = 0;
    for (int i = elem + 1; i < n; i++)
        summ += M[i][elem] * M[i][elem];
    return sqrt(summ);
}

double** findH(double* v1, double* v2, double denominator) {
    double** res = new double* [n];
    for (int i = 0; i < n; ++i) {
        res[i] = new double[n];
        for (int j = 0; j < n; ++j) {
            res[i][j] = - v1[i] * v2[j] * 2 / denominator;
        }
        res[i][i] = 1 + res[i][i];
    }
    return res;
}

double vector_multiple(double* v1, double* v2) {
    double res = 0;
    for (int i = 0; i < n; ++i) {
        res += v1[i] * v2[i];
    }
    return res;
}

double* get_vect(double** M, int col) {
    double* v = new double[n];
    double summ = 0;
    for (int i = 0; i < col; i++) {
        v[i] = 0;
    }
    for (int i = col; i < n; i++) {
        v[i] = M[i][col];
        summ += M[i][col] * M[i][col];
    }
    v[col] += sgn(v[col]) * sqrt(summ);
    return v;
}

int main() {
    ifstream inFile;

    inFile.open("matrix.txt");
    inFile >> eps;
    inFile >> n;
    double** A = new double* [n];
    double** Q = new double* [n];
    double** newH = new double* [n];
    double* v = new double[n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        Q[i] = new double[n];
        newH[i] = new double[n];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
        {
            Q[i][j] = 0;
            inFile >> A[i][j];
            cout << A[i][j] << " ";
        }
        Q[i][i] = 1;
        cout << "\n";
    }
    cout << "\n";

    double den, norm = 0, delta = -1, val, prev_val, under_root;
    complex* lambdas = new complex[n];
    //v = get_vect(A, 0);
    //den = vector_multiple(v, v);    //нахождение знаменателя
    //newH = findH(v, v, den); //нахождение H

    for (int k = 0; k < n; k++) {
        complex num, num2, old_num, old_num2;
        num.real = 0; num.im = 0;
        num2.real = 0; num2.im = 0;
        
        while (true) {
            for (int i = 0; i < n - 1; i++) {
                v = get_vect(A, i);
                den = vector_multiple(v, v);    //нахождение знаменателя
                newH = findH(v, v, den); //нахождение H
                A = multiple(newH, A);
                Q = multiple(Q, newH);
            }

            A = multiple(A, Q); //умножение A и Q
            ones(Q); // приведение Q к единичному виду
            print(A);

            if (find_norm(A, k) < eps) {
                lambdas[k].im = 0;
                lambdas[k].real = A[k][k];
                break;
            }

            old_num.real = num.real; old_num.im = num.im;
            old_num2.real = num2.real; old_num2.im = num2.im;
            under_root = A[k][k] * A[k][k] + A[k + 1][k + 1] * A[k + 1][k + 1] - 2 * A[k][k] * A[k + 1][k + 1] + 4 * A[k + 1][k] * A[k][k + 1];
            if (under_root < 0) {
                num.real = (A[k][k] + A[k + 1][k + 1]) / 2;
                num.im = sqrt(abs(under_root))/2;
                num2.real = (A[k][k] + A[k + 1][k + 1]) / 2;
                num2.im = -sqrt(abs(under_root))/2;
            }
            else {
                num.real = (A[k][k] + A[k + 1][k + 1]) / 2 + sqrt(under_root) / 2;;
                num.im = 0;
                num2.real = (A[k][k] + A[k + 1][k + 1]) / 2 - sqrt(under_root) / 2;
                num2.im = 0;
            }

            if (abs(num.real - old_num.real) < eps && abs(num2.real - old_num2.real) < eps) {
                lambdas[k] = num;
                lambdas[k + 1] = num2;
                k++;
                break;
            }
        }
    }

    for (int i = 0; i < n; i++) {
        cout << lambdas[i].real << " + " << lambdas[i].im << " * i || ";
    }
}
