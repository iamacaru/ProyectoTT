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
        friend Matrix operator*(const Matrix& matrix, double scalar);
        friend Matrix operator/(const Matrix& matrix, double scalar);
        friend Matrix operator+(const Matrix& matrix, double scalar);


        double& operator()(int i, int j) const;
 
        void print();

        double norm() const;
        Matrix cross(const Matrix& matrix2) const;

        int rows() const;
        int columns() const;

        Matrix traspuesta() const;
        Matrix inversa() const;
        Matrix subMatrix(int f, int l, int i) const;
        Matrix concatenar(const Matrix& matrix2) const;
        bool isEqual(const Matrix& matrix2, double TOL_) const;

        static Matrix identidad(int size);

 
    private:
        void initMatrix();
 
    private:
        int fil;
        int col;
        double **matrix;
};

#endif
