#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>

using namespace std;

int n;
double** A;

void swap_rows(double** M, int row1, int row2) {
    for (int i = 0; i < n; i++) {
        swap(M[row1][i], M[row2][i]);
    }
}

void print(double** M) {
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) {
            cout << M[i][j] << setw(15);
        }
        cout << "\n";
    }
    cout << "\n";
}

void pulDecomp(double** A, double** E, int* p, double* b) {
    double det = 1;
    int count = 0;
    double* z = new double[n], * x = new double[n];
    double** C = new double* [n];
    for (int i = 0; i < n; i++) {
        z[i] = 0;
        x[i] = 0;
        C[i] = new double[n];
    }

    
    for (int i = 0; i < n; i++) {
        //поиск опорного элемента
        double SuppElem = 0;
        int baseIndex = -1;
        for (int row = i; row < n; row++) {
            if (fabs(A[row][i]) > SuppElem) {
                SuppElem = fabs(A[row][i]);
                baseIndex = row;
            }
        }
        if (SuppElem != 0) {
            //меняем местами i-ю строку и строку с опорным элементом
            swap(p[baseIndex], p[i]);
            swap(b[baseIndex], b[i]);
            swap_rows(A, baseIndex, i);
            count++;
        }
    }
    
    print(A);

    for (int i = 1; i < n; i++)
        A[i][0] = A[i][0] / A[0][0];


    for (int main = 0; main < n-1; main++) {
        for (int k = main + 1; k < n; k++) {
            A[k][k] = A[k][k] - A[main][k] * A[k][main];
            b[k] = b[k] - b[main] * A[k][main];
        for (int i = 0; i < n; i++) {
            if (i > k) {
                A[k][i] = A[k][i] - A[k][main] * A[main][i];
                A[i][k] = A[i][k] - A[i][main] * A[main][k];

                if (k == main + 1) {
                    A[i][main + 1] = A[i][main + 1] / A[main + 1][main + 1];
                }
            } else
                E[k][i] = E[k][i] - A[k][main] * E[main][i];
        }
        }
        det *= A[main][main];
    }
    det *= A[n - 1][n - 1] * pow(-1, count);;
    cout << "det(A) = " << det << "\n";
    cout << "LU matrix:" << "\n";
    print(A);
    cout << "Reversed L matrix:" << "\n";
    print(E);

    for (int i = n - 1; i > -1; i--) {
        double ui = A[i][i], s = 0;
        for (int k = i + 1; k < n; k++)
            s += x[k] * A[i][k];
        x[i] = (b[i] - s) / ui;
    }

    cout << "X = ( ";
    for (int i = 0; i < n; i++) {
        cout << x[i] << " ";
    }
    cout << ")\n";
    
    for (int i = 0; i < n; i++) {
        for (int j = n - 1; j > -1; j--) {
            double ui = A[j][j];
            double s = 0;
            for (int k = j + 1; k < n; k++) {
                s += E[k][p[i]] * A[j][k];
            }
            E[j][p[i]] = (E[j][p[i]] - s) / ui;
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[j][p[i]] = E[j][i];
        }
    }
    cout << "reversed A matrix:" << "\n";
    print(C);
    
}

int main() {
    ifstream inFile;
    inFile.open("matrix.txt");
    inFile >> n;
    double** E = new double* [n];
    A = new double* [n];
    double* b = new double[n];
    int* p = new int[n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        E[i] = new double[n];
        p[i] = i;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
        {
            inFile >> A[i][j];
            E[i][j] = 0;
            cout << A[i][j] << " ";
        }
        E[i][i] = 1;
        cout << "\n";
    }

    for (int i = 0; i < n; i++) {
        inFile >> b[i];
        cout << b[i] << " ";
    }
    cout << "\n\n";

    pulDecomp(A, E, p, b);
}
