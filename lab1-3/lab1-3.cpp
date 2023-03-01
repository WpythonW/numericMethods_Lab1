#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

ifstream inFile;
int n;
double Cnorm = 0, CnormVect = 0, eps;

void input_matrix(double** A, double* b) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
        {
            inFile >> A[i][j];
            //cout << A[i][j] << " ";
        }
        //cout << "\n";
    }

    for (int i = 0; i < n; i++) {
        inFile >> b[i];
        //cout << b[i] << " ";
    }
    //cout << "\n";
}

void check_conditions(double** A, double* b) {
    double summ1 = 0, norm1 = 0, summ2 = 0, Csumm = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Csumm += abs(A[i][j]);
            summ2 += A[i][j] * A[i][j];
            summ1 += abs(A[j][i]);
        }
        if (norm1 < summ1) {
            norm1 = summ1;
        }
        if (Cnorm < Csumm) {
            Cnorm = Csumm;
        }
        summ1 = 0; Csumm = 0;
        if (CnormVect < abs(b[i]))
            CnormVect = abs(b[i]);
    }
    //cout << norm1 << " " << sqrt(summ2) << " " << Cnorm;

    if (norm1 >= 1 && sqrt(summ2) >= 1 && Cnorm >= 1)
    {
        cout << "norm > 1 => Sufficient condition not met\n";
        exit(0);
    }
}

void print(double* X) {
    for (int i = 0; i < n; i++) {
        cout << X[i] << " ";
    }
    cout << "\n";
}

double* iterations_method(double** A, double* b) {
    int k = 0;
    double* X = new double[n], * oldX = new double[n], norm = 0;
    copy(b, b + n, oldX);
    double epsK, s = 0;
    do {
        k++;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                s += A[i][j] * oldX[j];
            }
            X[i] = s + b[i];
            s = 0;

            if (norm < abs(X[i] - oldX[i]))
                norm = abs(X[i] - oldX[i]);
        }
        epsK = (Cnorm / (1 - Cnorm)) * norm;
        copy(X, X + n, oldX);
        norm = 0;
        //cout << epsK << "\n";
    } while (epsK > eps);

    return X;
}

double* seidel_method(double** A, double* b) {
    int k = 0;
    double* X = new double[n], * oldX = new double[n], norm = 0, norm_upper = 0, summ = 0;
    copy(b, b + n, oldX);
    double epsK, s = 0;
    do {
        k++;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j < i)
                    s += A[i][j] * X[j];
                else {
                    s += A[i][j] * oldX[j];
                    if (norm_upper == 0) summ += abs(A[i][j]);
                }

            }
            X[i] = s + b[i];
            s = 0;

            if (norm < abs(X[i] - oldX[i]))
                norm = abs(X[i] - oldX[i]);
            if (norm_upper < summ && norm_upper == 0)
                norm_upper = summ;
        }
        epsK = (norm_upper / (1 - Cnorm)) * norm;
        copy(X, X + n, oldX);
        norm = 0;
        //cout << epsK << "\n";
        //print(X);
    } while (epsK > eps);

    return X;
}

int main() {
    inFile.open("matrix.txt");
    inFile >> eps;
    inFile >> n;
    double** A = new double* [n];
    double k;
    for (int i = 0; i < n; i++)
        A[i] = new double[n];

    double* b = new double[n], * X_iter = new double[n], * X_Seidel = new double[n];

    input_matrix(A, b); //Ввод матрицы

    for (int i = 0; i < n; i++) {
        double elem = A[i][i];
        for (int j = 0; j < n; j++) {
            A[i][j] = -A[i][j] / elem;
            if (i == j)
                A[i][j] = 0;
            //cout << A[i][j] << " ";
        }
        b[i] = b[i] / elem;
    }
    // приведение к эквивалентному виду
    check_conditions(A, b); //проверка на достаточное условие

    k = (log(eps) - log(CnormVect) + log(1 - Cnorm)) / log(Cnorm) - 1;
    cout << "a priori estimate of the number of iterations k: " << k << "\n";

    X_iter = iterations_method(A, b);
    cout << "Vector X found by iterations method\n";
    print(X_iter);

    cout << "\n";

    X_Seidel = seidel_method(A, b);
    cout << "Vector X found by Seidel method\n";
    print(X_Seidel);
}
