#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>

using namespace std;

int n;

void swap_rows(double** M, int row1, int row2) {
    for (int i = 0; i < n; i++) {
        swap(M[row1][i], M[row2][i]);
    }
}

void print(double** M) {
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) {
            cout << round(M[i][j]*10000)/10000 << setw(10);
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

    for (int i = 0; i < n; i++){
        double elem = A[i][i];
        double coeff;
        for (int j = i + 1; j < n; j++) {
            coeff = A[j][i] / elem;
            for (int k = 0; k < n; k++) {
                if (k >= i)
                    A[j][k] = A[j][k] - A[i][k] * coeff; //Условие здесь, чтобы не залезать на матрицу L
                E[j][k] = E[j][k] - E[i][k] * coeff;
            }
            A[j][i] = coeff;
        }
        det *= A[i][i];
    }
    det *= pow(-1, count);

    //print(A);

    for (int i = 0; i < n; i++) {
        double s = 0;
        for (int k = 0; k < i; k++)
            s += z[k] * A[i][k];
        z[i] = b[i] - s;
        //cout << z[i] << " ";
    }

    for (int i = n - 1; i > -1; i--) {
        double ui = A[i][i], s = 0;
        for (int k = i + 1; k < n; k++)
            s += x[k] * A[i][k];
        x[i] = (z[i] - s) / ui;
    }

    cout << "X = ( ";
    for (int i = 0; i < n; i++) {
        cout << x[i] << " ";
    }
    cout << ")\n";
    
    for (int i = 0; i < n; i++) {
        for (int j = n-1; j > -1; j--) {
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
    print(C);
}

int main() {
    ifstream inFile;
    inFile.open("matrix.txt");
    inFile >> n;
    double** A = new double* [n], **E = new double* [n];
    double* b = new double[n];
    int *p = new int[n];
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
