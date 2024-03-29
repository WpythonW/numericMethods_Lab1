#include <iostream>
#include <fstream>

using namespace std;

int main() {
    int n;
    bool cond = false;
    ifstream inFile;
    inFile.open("matrix.txt");
    inFile >> n;
    double** A = new double* [n];
    for (int i = 0; i < n; i++)
        A[i] = new double[n];
    double* d = new double[n], * P = new double[n], * Q = new double[n], * X = new double[n];
    char ch;
    for (int i = 0; i < n; i++) {
        inFile >> A[i][0];
        inFile >> A[i][1];
        inFile >> A[i][2];
        if ((abs(A[i][0]) + abs(A[i][2]) > abs(A[i][1])) || ((A[i][0] == 0 || A[i][2] == 0) && i != 0 && i != n - 1)) {
            cond = true;
        }
        cout << A[i][0] << " " << A[i][1] << " " << A[i][2] << "\n";
    }

    for (int i = 0; i < n; i++) {
        inFile >> d[i];
        cout << d[i] << " ";
    }
    cout << "\n";

    P[0] = -A[0][2] / A[0][1];
    Q[0] = d[0] / A[0][1];

    for (int i = 1; i < n; i++)
    {
        P[i] = -A[i][2] / (A[i][1] + A[i][0] * P[i - 1]);
        Q[i] = (d[i] - A[i][0] * Q[i - 1]) / (A[i][1] + A[i][0] * P[i - 1]);

    }

    X[n - 1] = Q[n - 1];

    for (int i = n - 2; i >= 0; i--)
    {
        X[i] = P[i] * X[i + 1] + Q[i];
    }
    cout << "\n";
    for (int i = 0; i < n; i++) {
        cout << X[i] << " ";
    }

    if (cond) cout << "\nSufficient condition not met";
}
