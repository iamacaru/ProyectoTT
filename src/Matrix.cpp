#include "../include/Matrix.h"

#include <iostream>
#include <iomanip>
#include <cmath>

Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}
 
Matrix::Matrix(int fil, int col, double v[], int n): fil(fil), col(col)
{
    initMatrix();
 
    int k = 0;
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}
 
Matrix::Matrix(const Matrix& m)
{
    *this = m;
}
 
Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];
 
    delete[] matrix;
}
 
void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];
 
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}
 
Matrix& Matrix::operator=(const Matrix& matrix2)
{
    this->fil = matrix2.rows();
    this->col = matrix2.columns();

    initMatrix();

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            this->matrix[i][j] = matrix2.matrix[i][j];

    return *this;
}
 
Matrix Matrix::operator+(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator-(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator*(const Matrix& matrix2)
{
    Matrix result(this->fil, matrix2.col);
 
    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }
 
    return result;
}

Matrix operator*(double scalar, const Matrix& matrix) {
    Matrix result(matrix.rows(), matrix.columns());
    for (int i = 0; i < matrix.rows(); i++) {
        for (int j = 0; j < matrix.columns(); j++) {
            result(i + 1, j + 1) = scalar * matrix(i + 1, j + 1);
        }
    }
    return result;
}

Matrix operator*(const Matrix& matrix, double scalar) {
    Matrix result(matrix.rows(), matrix.columns());
    for (int i = 0; i < matrix.rows(); i++) {
        for (int j = 0; j < matrix.columns(); j++) {
            result(i + 1, j + 1) = scalar * matrix(i + 1, j + 1);
        }
    }
    return result;
}

Matrix operator/(const Matrix& matrix, double scalar) {
    if (scalar == 0.0) {
        throw std::invalid_argument("Division by zero is not allowed.");
    }

    Matrix result(matrix.rows(), matrix.columns());

    for (int i = 0; i < matrix.rows(); i++) {
        for (int j = 0; j < matrix.columns(); j++) {
            result(i + 1, j + 1) = matrix(i + 1, j + 1) / scalar;
        }
    }

    return result;
}

Matrix operator+(const Matrix& matrix, double scalar) {
    Matrix result(matrix.rows(), matrix.columns());

    for (int i = 0; i < matrix.rows(); i++) {
        for (int j = 0; j < matrix.columns(); j++) {
            result(i + 1, j + 1) = matrix(i + 1, j + 1) + scalar;
        }
    }

    return result;
}
 
double& Matrix::operator()(const int i, const int j) const
{
    return matrix[i-1][j-1];
}
 
void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

double Matrix::norm() const{
    double sum = 0.0;

    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            sum += matrix[i][j] * matrix[i][j];
        }
    }

    return sqrt(sum);
}

Matrix Matrix::cross(const Matrix& matrix2) const {
    double x = matrix[0][1] * matrix2(1,3) - matrix[0][2] * matrix2(1,2);
    double y = matrix[0][2] * matrix2(1,1) - matrix[0][0] * matrix2(1,3);
    double z = matrix[0][0] * matrix2(1,2) - matrix[0][1] * matrix2(1,1);

    Matrix result(1, 3);
    result(1, 1) = x;
    result(1, 2) = y;
    result(1, 3) = z;

    return result;
}

int Matrix::rows() const {
    return fil;
}

int Matrix::columns() const {
    return col;
}

Matrix Matrix::traspuesta() const {
    double values[fil*col];

    int k = 0;

    for (int i = 0; i < col; i++) {
        for (int j = 0; j < fil; j++) {
            values[k] = matrix[j][i];
            k++;
        }
    }

    return {col, fil, values, col * fil};
}

Matrix Matrix::inversa() const {
    int n = fil;
    Matrix A(*this);
    Matrix I = Matrix::identidad(n);

    // Gaussian elimination
    for (int i = 1; i <= n; i++) {
        // Search for maximum in this column
        double maxEl = std::abs(A(i, i));
        int maxRow = i;
        for (int k = i + 1; k <= n; k++) {
            if (std::abs(A(k, i)) > maxEl) {
                maxEl = std::abs(A(k, i));
                maxRow = k;
            }
        }

        // Swap maximum row with current row (pivot)
        for (int k = 1; k <= n; k++) {
            std::swap(A(maxRow, k), A(i, k));
            std::swap(I(maxRow, k), I(i, k));
        }

        // Make all rows below this one 0 in the current column
        for (int k = i + 1; k <= n; k++) {
            double c = -A(k, i) / A(i, i);
            for (int j = 1; j <= n; j++) {
                if (i == j) {
                    A(k, j) = 0;
                } else {
                    A(k, j) += c * A(i, j);
                }
                I(k, j) += c * I(i, j);
            }
        }
    }

    // Solve equation Ax = I for an upper triangular matrix A
    for (int i = n; i >= 1; i--) {
        for (int j = 1; j <= n; j++) {
            I(i, j) /= A(i, i);
        }
        for (int k = i - 1; k >= 1; k--) {
            double c = -A(k, i);
            for (int j = 1; j <= n; j++) {
                I(k, j) += c * I(i, j);
            }
        }
    }

    return I;
}



Matrix Matrix::subMatrix(int f, int l, int i) const{
    Matrix r(1, l - f + 1);

    int k = 1;
    for (int j = f; j <= l; j++) {
        r(1, k) = matrix[i - 1][j - 1];
        k++;
    }

    return r;
}

Matrix Matrix::concatenar(const Matrix& matrix2) const {
    Matrix r(1, this->col + matrix2.col);

    for (int j = 0; j < this->col; ++j) {
        r(1, j + 1) = matrix[0][j];
    }

    for (int j = 1; j <= matrix2.col; ++j) {
        r(1, this->col + j) = matrix2(1, j);
    }

    return r;
}

bool Matrix::isEqual(const Matrix& matrix2, double TOL_) const {
    if (this->fil != matrix2.fil || this->col != matrix2.col) {
        return false;
    }

    for (int i = 0; i < this->fil; ++i) {
        for (int j = 0; j < this->col; ++j) {
            if (fabs(this->matrix[i][j] - matrix2.matrix[i][j]) > TOL_) {
                return false;
            }
        }
    }

    return true;
}

Matrix Matrix::identidad(int size) {
    Matrix result(size, size);
    for (int i = 1; i <= size; ++i) {
        for (int j = 1; j <= size; ++j) {
            result(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }
    return result;
}