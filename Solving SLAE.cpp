#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

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
            cout << fixed << setprecision(2) << v[i] << endl;
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

    Matrix operator + (const Matrix& b) {
        Matrix res(m, b.n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                res.v[i][j] = v[i][j] + b.v[i][j];
            }
        }
        return res;
    }

    Matrix operator - (const Matrix& b) {
        Matrix res(m, b.n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                res.v[i][j] = v[i][j] - b.v[i][j];
            }
        }
        return res;
    }

    Matrix operator * (const Matrix& b) {
        Matrix res(m, b.n);
        for (int i = 0; i < res.m; i++) {
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
                cout << fixed << setprecision(2) << v[i][j] << " ";
            }
            cout << endl;
        }
    }
};
void outputIdentityMatrix(int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i == j)
                cout << "1.00 ";
            else
                cout << "0.00 ";
        }
        cout << endl;
    }
}
class SquareMatrix : public Matrix {
public:
    SquareMatrix(int a) : Matrix(a, a) {}

    double determinant() {
        if (m != n) {
            cout << "Error: Determinant is only defined for square matrices." << endl;
            return 0.0;
        }

        SquareMatrix temp(*this); // Make a copy of the matrix to avoid modifying the original

        double det = 1.0;

        for (int i = 0; i < n; ++i) {
            if (temp.v[i][i] == 0) {
                // Find a row with a nonzero pivot element to swap
                bool found = false;
                for (int k = i + 1; k < n; ++k) {
                    if (temp.v[k][i] != 0) {
                        swap(temp.v[i], temp.v[k]);
                        det = -det; // Swap changes the sign of the determinant
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    // If no such row is found, determinant is zero
                    return 0.0;
                }
            }

            // Perform row operations to make elements below the pivot zero
            for (int j = i + 1; j < n; ++j) {
                double factor = temp.v[j][i] / temp.v[i][i];
                for (int k = i; k < n; ++k) {
                    temp.v[j][k] -= factor * temp.v[i][k];
                }
            }

            // Multiply the determinant by the pivot element
            det *= temp.v[i][i];
        }

        return det;
    }
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix(int a) : SquareMatrix(a) {
        for (int i = 0; i < a; i++) {
            v[i][i] = 1;
        }
    }
};

class EliminationMatrix : public IdentityMatrix {
public:
    EliminationMatrix(int a, int x, int y, Matrix matrix) : IdentityMatrix(a) {
        if (matrix.v[y - 1][y - 1] != 0) {
            v[x - 1][y - 1] = -(double(matrix.v[x - 1][y - 1] / matrix.v[y - 1][y - 1]));
        }
    }
};

class PermutationMatrix : public IdentityMatrix {
public:
    PermutationMatrix(int a, int r1, int r2) : IdentityMatrix(a) {
        v[r1 - 1][r1 - 1] = 0;
        v[r1 - 1][r2 - 1] = 1;
        v[r2 - 1][r2 - 1] = 0;
        v[r2 - 1][r1 - 1] = 1;
    }
};

void inverse(SquareMatrix* A, ColumnVector* b) {
    int n = A->n;
    int step = 1;
    Matrix* I = new IdentityMatrix(n);
    for (int i = 1; i < n; i++) {
        int maxPivot = i - 1;
        for (int j = i - 1; j < n; j++) {
            if (fabs(A->v[j][i - 1]) > fabs(A->v[maxPivot][i - 1])) {
                maxPivot = j;
            }
        }
        if (fabs(A->v[maxPivot][i - 1]) < 1e-9) {
            cout << "Error: matrix A is singular" << endl;
            return;
        }
        if (maxPivot + 1 != i) {
            Matrix* P = new PermutationMatrix(n, maxPivot + 1, i);
            *(Matrix*)A = *P * *A;
            (*I) = *P * *I;
            *b = *P * *b; // Apply permutation to vector b
            // Output permutation step
            cout << "step #" << step << ": permutation" << endl;
            A->output();
            b->output();
            step++;
        }
        for (int j = i + 1; j <= n; j++) {
            if (A->v[j - 1][i - 1] != 0) {
                Matrix* E = new EliminationMatrix(n, j, i, *((Matrix*)A));
                *(Matrix*)A = *E * *A;
                (*I) = *E * *I;
                *b = *E * *b; // Apply elimination to vector b
                // Output elimination step
                cout << "step #" << step << ": elimination" << endl;
                A->output();
                b->output();
                step++;
            }
        }
    }
    for (int i = n; i > 1; i--) {
        for (int j = i - 1; j >= 1; j--) {
            if (A->v[j - 1][i - 1] != 0) {
                double div = -(A->v[j - 1][i - 1] / A->v[i - 1][i - 1]);
                Matrix* temp = new IdentityMatrix(n);
                temp->v[j - 1][i - 1] = div;
                *(Matrix*)A = *temp * *A;
                (*I) = *temp * *I;
                *b = *temp * *b; // Apply elimination to vector b
                // Output elimination step
                cout << "step #" << step << ": elimination" << endl;
                A->output();
                b->output();
                step++;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        Matrix* temp = new IdentityMatrix(n);
        if (A->v[i][i] != 0) {
            temp->v[i][i] = (1 / A->v[i][i]);
            *(Matrix*)A = *temp * *A;
            (*I) = *temp * *I;
        }
    }
    *(Matrix*)A = *I;
}
void inputVector(ColumnVector& v) {
    for (int i = 0; i < v.n; ++i) {
        cin >> v.v[i];
    }
}

int main() {
    int size;
    cin >> size;

    SquareMatrix A(size);
    ColumnVector b(size);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cin >> A.v[i][j];
        }
    }

    cin >> size;
    inputVector(b);

    // Check determinant
    double det = A.determinant();
    if (det != 0.0) {
        // Compute inverse
        cout<<"Gaussian process:"<<endl;
        ColumnVector temp = b;
        inverse(&A,&temp);
        // Calculate solution
        ColumnVector x = A * b;

        cout << "Diagonal normalization:" << endl;
        outputIdentityMatrix(size);
        x.output();
        // Output solution
        cout << "Result:" << endl;
        x.output();
    } else {
        cout << "Error: matrix A is singular" << endl;
    }

    return 0;
}