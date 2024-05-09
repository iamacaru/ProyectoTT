#ifndef _MATRIX_
#define _MATRIX_

class Matrix
{
    public:
        Matrix(int fil, int col);
        Matrix(int fil, int col, double v[], int n);
        Matrix(const Matrix& m);
        ~Matrix();
 
        Matrix& operator=(const Matrix& matrix2);
        Matrix  operator+(const Matrix& matrix2);
        Matrix  operator-(const Matrix& matrix2);
        Matrix  operator*(const Matrix& matrix2);
        friend Matrix operator*(double scalar, const Matrix& matrix);
        friend Matrix operator/(const Matrix& matrix, double scalar);


        double& operator()(int i, int j) const;
 
        void print();

        double norm() const;
        double dot(const Matrix& matrix2) const;

        int rows() const;
        int columns() const;


 
    private:
        void initMatrix();
 
    private:
        int fil;
        int col;
        double **matrix;
};

#endif
