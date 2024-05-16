#include <bits/stdc++.h>

using namespace std;

class MatrixBase {
public:
    virtual void input() = 0;
    virtual void output() const = 0;
    virtual MatrixBase* add(const MatrixBase& other) const = 0;
    virtual MatrixBase* subtract(const MatrixBase& other) const = 0;
    virtual MatrixBase* multiply(const MatrixBase& other) const = 0;
    virtual MatrixBase* transpose() const = 0;
    virtual ~MatrixBase() {}
};

class Matrix : public MatrixBase {
private:
    vector<vector<int>> data;
    int size;

public:
    Matrix(int n) : size(n) {
        data.resize(size, vector<int>(size, 0));
    }

    void input() override {
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                cin >> data[i][j];
            }
        }
    }

    void output() const override {
        for (const auto& row : data) {
            for (const int& element : row) {
                cout << element << " ";
            }
            cout << endl;
        }
    }

    int getSize() const {
        return size;
    }

    MatrixBase* add(const MatrixBase& other) const override {
        const Matrix& otherMatrix = dynamic_cast<const Matrix&>(other);
        if (size != otherMatrix.size) {
            cout << "Error: the dimensional problem occurred" << endl;
            return nullptr;
        }
        Matrix* result = new Matrix(size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                result->data[i][j] = data[i][j] + otherMatrix.data[i][j];
            }
        }
        return result;
    }

    MatrixBase* subtract(const MatrixBase& other) const override {
        const Matrix& otherMatrix = dynamic_cast<const Matrix&>(other);
        if (size != otherMatrix.size) {
            cout << "Error: the dimensional problem occurred" << endl;
            return nullptr;
        }
        Matrix* result = new Matrix(size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                result->data[i][j] = data[i][j] - otherMatrix.data[i][j];
            }
        }
        return result;
    }

    MatrixBase* multiply(const MatrixBase& other) const override {
        const Matrix& otherMatrix = dynamic_cast<const Matrix&>(other);
        if (size != otherMatrix.size) {
            cout << "Error: the dimensional problem occurred" << endl;
            return nullptr;
        }
        Matrix* result = new Matrix(size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                for (int k = 0; k < size; ++k) {
                    result->data[i][j] += data[i][k] * otherMatrix.data[k][j];
                }
            }
        }
        return result;
    }

    MatrixBase* transpose() const override {
        Matrix* result = new Matrix(size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                result->data[i][j] = data[j][i];
            }
        }
        return result;
    }
};

istream& operator>>(istream& in, MatrixBase& matrix) {
    matrix.input();
    return in;
}

ostream& operator<<(ostream& out, const MatrixBase& matrix) {
    matrix.output();
    return out;
}

int main() {
    int size;

    // Input matrix A
    cin >> size;
    Matrix A(size);
    cin >> A;

    // Input matrix B
    cin >> size;
    Matrix B(size);
    cin >> B;

    // Input matrix C
    cin >> size;
    Matrix C(size);
    cin >> C;

    // Perform operations and output results or error messages

    // Matrix Addition
    MatrixBase* D = A.add(B);
    if (D != nullptr) {
        cout << *D;
        delete D;
    }

    // Matrix Subtraction
    MatrixBase* E = B.subtract(A);
    if (E != nullptr) {
        cout << *E;
        delete E;
    }

    // Matrix Multiplication
    MatrixBase* F = C.multiply(A);
    if (F != nullptr) {
        cout << *F;
        delete F;
    }

    // Matrix Transpose
    MatrixBase* G = A.transpose();
    if (G != nullptr) {
        cout << *G;
        delete G;
    }

    return 0;
}