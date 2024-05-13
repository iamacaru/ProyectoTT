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
    Matrix result(fil, col);
 
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

double Matrix::dot(const Matrix& matrix2) const {
    if (fil != 1 || matrix2.rows() != 1) {
        throw std::invalid_argument("Both matrices must have one row.");
    }

    if (col != matrix2.columns()) {
        throw std::invalid_argument("Matrices must have the same number of columns.");
    }

    double result = 0.0;

    for (int j = 0; j < col; j++) {
        result += matrix[0][j] * matrix2.matrix[0][j];
    }

    return result;
}

Matrix Matrix::cross(const Matrix& matrix2) const {
    double x = matrix[0][1] * matrix2(0,2) - matrix[0][2] * matrix2(0,1);
    double y = matrix[0][2] * matrix2(0,0) - matrix[0][0] * matrix2(0,2);
    double z = matrix[0][0] * matrix2(0,1) - matrix[0][1] * matrix2(0,0);

    Matrix result(1, 3);
    result(0, 0) = x;
    result(0, 1) = y;
    result(0, 2) = z;

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

Matrix Matrix::subMatrix(int f, int l, int w, bool b) const{
    int k = 1;
    if (b) {
        Matrix r(1, l - f + 1);
        for (int i = f; i <= l; i++) {
            r(1, k) = matrix[w - 1][i - 1];
            k++;
        }

        return r;
    } else {
        Matrix r(l - f + 1, 1);
        for (int i = f; i <= l; i++) {
            r(k, 1) = matrix[i - 1][w - 1];
            k++;
        }

        return r;
    }
}