#include <bits/stdc++.h>
using namespace std;

class ColumnVector {
public:
    int n;
    vector<double> v;
    ColumnVector(int m) {
        v = vector<double>(m, 0);
        n = m;
    }

    void operator = (ColumnVector a) {
        for (int i = 0; i < n; i++) {
            v[i] = a.v[i];
        }
    }

    ColumnVector operator + (const ColumnVector& b) {
        ColumnVector res(n);
        for (int i = 0; i < n; i++) {
            res.v[i] = v[i] + b.v[i];
        }
        return res;
    }

    ColumnVector operator - (const ColumnVector& b) {
        ColumnVector res(n);
        for (int i = 0; i < n; i++) {
            res.v[i] = v[i] - b.v[i];
        }
        return res;
    }

    void output() {
        for (int i = 0; i < n; i++) {
            cout << fixed << setprecision(4) << v[i] << endl;
        }
    }
};

class Matrix {
public:
    int m, n;
    vector<vector<double>> v;
    Matrix(int a, int b) {
        m = a;
        n = b;
        v = vector<vector<double>>(m, vector<double>(n, 0));
    }

    void operator = (Matrix a) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                v[i][j] = a.v[i][j];
            }
        }
    }

    Matrix operator * (const Matrix& b) {
        Matrix res(m, b.n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < b.n; j++) {
                double sum = 0;
                for (int k = 0; k < n; k++) {
                    sum += (v[i][k] * b.v[k][j]);
                }
                if (fabs(sum) <= 1e-4) sum = 0;
                res.v[i][j] = sum;
            }
        }
        return res;
    }

    ColumnVector operator * (const ColumnVector& b) {
        ColumnVector res(m);
        for (int i = 0; i < m; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                sum += (v[i][j] * b.v[j]);
            }
            if (fabs(sum) <= 1e-4) sum = 0;
            res.v[i] = sum;
        }
        return res;
    }

    Matrix T() {
        Matrix res(n, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                res.v[j][i] = v[i][j];
            }
        }
        return res;
    }

    void output() {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                cout << fixed << setprecision(4) << v[i][j] << " ";
            }
            cout << endl;
        }
    }
};

void inverse(Matrix* A) {
    int n = A->n;
    Matrix* I = new Matrix(n, n);
    for (int i = 0; i < n; i++) {
        (*I).v[i][i] = 1;
    }
    for (int i = 0; i < n; i++) {
        int maxPivot = i;
        for (int j = i + 1; j < n; j++) {
            if (fabs(A->v[j][i]) > fabs(A->v[maxPivot][i])) {
                maxPivot = j;
            }
        }
        if (fabs(A->v[maxPivot][i]) < 1e-9) {
            cout << "Error: matrix A is singular" << endl;
            return;
        }
        if (maxPivot != i) {
            swap(A->v[i], A->v[maxPivot]);
            swap((*I).v[i], (*I).v[maxPivot]);
        }
        for (int j = i + 1; j < n; j++) {
            double factor = A->v[j][i] / A->v[i][i];
            for (int k = 0; k < n; k++) {
                A->v[j][k] -= factor * A->v[i][k];
                (*I).v[j][k] -= factor * (*I).v[i][k];
            }
        }
    }
    for (int i = n - 1; i > 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            double factor = A->v[j][i] / A->v[i][i];
            for (int k = 0; k < n; k++) {
                A->v[j][k] -= factor * A->v[i][k];
                (*I).v[j][k] -= factor * (*I).v[i][k];
            }
        }
    }
    for (int i = 0; i < n; i++) {
        double div = 1 / A->v[i][i];
        for (int j = 0; j < n; j++) {
            A->v[i][j] *= div;
            (*I).v[i][j] *= div;
        }
    }
    *A = *I;
}

ColumnVector LeastSquareApproximation(int n, int deg, vector<double> ts, vector<double> bs) {
    Matrix* A = new Matrix(n, deg + 1);
    ColumnVector B(n);
    for (int i = 0; i < deg + 1; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0)
                A->v[j][i] = 1;
            else if (i == 1)
                A->v[j][i] = ts[j];
            else
                A->v[j][i] = pow(ts[j], i);
        }
    }
    for (int i = 0; i < n; i++) {
        B.v[i] = bs[i];
    }
    Matrix* ATA = new Matrix(deg + 1, deg + 1);
    *ATA = A->T() * *A;
    cout << "A:" << endl;
    A->output();
    cout << "A_T*A:" << endl;
    ATA->output();
    Matrix* Ainv = new Matrix(deg + 1, deg + 1);
    *Ainv = *ATA;
    inverse(Ainv);
    ColumnVector ATb(deg + 1);
    ATb = A->T() * B;
    cout << "(A_T*A)_-1:" << endl;
    Ainv->output();
    ColumnVector x(deg + 1);
    x = (*Ainv) * ATb;
    ColumnVector result = (*Ainv) * ATb;
    cout << "A_T*b:" << endl;
    ATb.output();
    cout << "x~:" << endl;
    result.output();
    return x;
}


int main() {
    int m, n;
    cin >> m;
    vector<double> ts(m);
    vector<double> bs(m);
    for (int i = 0; i < m; ++i) {
        cin >> ts[i] >> bs[i];
    }
    cin >> n;
    ColumnVector x = LeastSquareApproximation(m, n, ts, bs);
    return 0;
}