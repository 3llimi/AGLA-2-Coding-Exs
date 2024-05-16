#include <bits/stdc++.h>

using namespace std;

class Matrix {
private:
    vector<vector<int>> data;
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

    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            cout << "Error: the dimensional problem occurred" << endl;
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
            cout << "Error: the dimensional problem occurred" << endl;
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
            cout << "Error: the dimensional problem occurred" << endl;
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

    Matrix& operator=(const Matrix& other) {
        if (this == &other) return *this;
        rows = other.rows;
        cols = other.cols;
        data = other.data;
        return *this;
    }
};

int main() {
    int rows, cols;
    // Input matrix A
    cin >> rows >> cols;
    Matrix A(rows, cols);
    A.input();
    // Input matrix B
    cin >> rows >> cols;
    Matrix B(rows, cols);
    B.input();
    // Input matrix C
    cin >> rows >> cols;
    Matrix C(rows, cols);
    C.input();

    // Perform operations and output results
    Matrix D = A + B;
    Matrix E = B - A;
    Matrix F = C * A;
    Matrix G = A.transpose();

    D.output();
    E.output();
    F.output();
    G.output();

    return 0;
}
