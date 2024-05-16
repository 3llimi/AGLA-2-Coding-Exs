#include <bits/stdc++.h>
using namespace std;

using namespace std;

class Matrix {
protected:
    int rows;
    int cols;

public:
    vector<vector<double>> data;

    Matrix(int r, int c) : rows(r), cols(c) {
        data.resize(rows, vector<double>(cols, 0));
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
            for (const double& element : row) {
                cout << fixed << setprecision(2) << element << " ";
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

        if (originalMatrix.data[y - 1][y - 1] == 0) {
            return;
        }

        double multiplier = -originalMatrix.data[x - 1][y - 1] / originalMatrix.data[y - 1][y - 1];

        data[x - 1][y - 1] = multiplier;
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
    cin >> size;

    Matrix A(size, size);
    A.input();

    Matrix temp = A;

    double determinant = 1.0;
    int stepCount = 0;
    int permCount = 0;

    for (int step = 0; step < size; ++step) {
        double maxVal = abs(temp.data[step][step]);
        int maxIndex = step;
        for (int i = step + 1; i < size; ++i) {
            if (abs(temp.data[i][step]) > maxVal) {
                maxVal = abs(temp.data[i][step]);
                maxIndex = i;
            }
        }

        if (maxVal == 0) {
            determinant = 0;
            break;
        }

        if (maxIndex != step) {
            PermutationMatrix P(size, step, maxIndex);
            temp = P * temp;
            cout << "step #" << ++stepCount << ": permutation" << endl;
            permCount++;
            temp.output();
        }

        for (int i = step + 1; i < size; ++i) {
            EliminationMatrix E(size, i + 1, step + 1, temp);
            temp = E * temp;
            cout << "step #" << ++stepCount << ": elimination" << endl;
            temp.output();
        }

        determinant *= temp.data[step][step];
    }
    if (permCount%2 != 0){
        determinant = -determinant;
    }
    cout << "result:" << endl;
    cout << fixed << setprecision(2) << determinant << endl;

    return 0;
}