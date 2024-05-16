#include <bits/stdc++.h>
using namespace std;

class Matrix {
protected:
    int rows;
    int cols;

public:
    Matrix(int r, int c) : rows(r), cols(c) {
        data.resize(rows, vector<int>(cols, 0));
    }

    void input() {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                cin >> data[i][j];
            }
        }
    }

    void output() const {
        for (const auto& row : data) {
            for (const int& element : row) {
                cout << element << " ";
            }
            cout << endl;
        }
    }

    int getRows() const {
        return rows;
    }

    int getCols() const {
        return cols;
    }

    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            return Matrix(0, 0);
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            return Matrix(0, 0);
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            return Matrix(0, 0);
        }
        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                for (int k = 0; k < cols; ++k) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    Matrix transpose() const {
        Matrix result(cols, rows);
        for (int i = 0; i < cols; ++i) {
            for (int j = 0; j < rows; ++j) {
                result.data[i][j] = data[j][i];
            }
        }
        return result;
    }
    vector<vector<int>> data;
};

class IdentityMatrix : public Matrix {
public:
    IdentityMatrix(int n) : Matrix(n, n) {
        for (int i = 0; i < n; ++i) {
            data[i][i] = 1;
        }
    }
};
class EliminationMatrix : public Matrix {
public:
    EliminationMatrix(int n, int x, int y, const Matrix& originalMatrix) : Matrix(n, n) {
        for (int i = 0; i < n; ++i) {
            data[i][i] = 1;
        }

        if (originalMatrix.data[y-1][y-1] == 0) {
            return;
        }

        int multiplier = -originalMatrix.data[x-1][y-1] / originalMatrix.data[y-1][y-1];

        data[x-1][y-1] = multiplier;
    }
};

class PermutationMatrix : public Matrix {
public:
    PermutationMatrix(int n, int i1, int i2) : Matrix(n, n) {
        for (int i = 0; i < n; ++i) {
            if (i == i1) {
                data[i][i1] = 0;
                data[i][i2] = 1;
            } else if (i == i2) {
                data[i][i2] = 0;
                data[i][i1] = 1;
            } else {
                data[i][i] = 1;
            }
        }
    }
};

int main() {
    int size;

    // Input matrix A
    cin >> size;
    Matrix A(size, size);
    A.input();

    // Output identity matrix of size 3x3
    IdentityMatrix I(3);
    I.output();

    // Output elimination matrix E21
    EliminationMatrix E(size, 2, 1, A);
    E.output();

    // Matrix B = E21 * A
    Matrix B = E * A;
    B.output();

    // Output permutation matrix P21
    PermutationMatrix P(size, 1, 0);
    P.output();

    // Matrix C = P21 * A
    Matrix C = P * A;
    C.output();

    return 0;
}