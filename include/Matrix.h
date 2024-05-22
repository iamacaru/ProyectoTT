#ifndef _MATRIX_
#define _MATRIX_

/*!
 * @file Matrix.h
 * @class Matrix
 * @brief A class for handling matrices and their operations
 */
class Matrix
{
    public:
         /*!
         * @brief Constructor that initializes a matrix with given dimensions
         *
         * @param fil Number of rows
         * @param col Number of columns
         */
        Matrix(int fil, int col);

        /*!
         * @brief Constructor that initializes a matrix with given dimensions and values from an array
         *
         * @param fil Number of rows
         * @param col Number of columns
         * @param v Array containing values to initialize the matrix
         * @param n Length of the array 'v'
         */
        Matrix(int fil, int col, double v[], int n);

        /*!
         * @brief Copy constructor
         *
         * @param m Another matrix to copy from
         */
        Matrix(const Matrix& m);

        //! Destructor.
        ~Matrix();

        /*!
         * @brief Assignment operator
         *
         * @param matrix2 Another matrix to assign from
         * @return Reference to the assigned matrix
         */
        Matrix& operator=(const Matrix& matrix2);

        /*!
         * @brief Addition operator
         *
         * @param matrix2 Another matrix to add
         * @return Resultant matrix after addition
         */
        Matrix  operator+(const Matrix& matrix2);

        /*!
         * @brief Subtraction operator
         *
         * @param matrix2 Another matrix to subtract
         * @return Resultant matrix after subtraction
         */
        Matrix  operator-(const Matrix& matrix2);

        /*!
         * @brief Multiplication operator
         *
         * @param matrix2 Another matrix to multiply
         * @return Resultant matrix after multiplication
         */
        Matrix  operator*(const Matrix& matrix2);

        /*!
         * @brief Scalar multiplication operator
         *
         * @param scalar Scalar value to multiply
         * @param matrix Matrix to multiply with the scalar
         * @return Resultant matrix after scalar multiplication
         */
        friend Matrix operator*(double scalar, const Matrix& matrix);

        /*!
         * @brief Scalar multiplication operator
         *
         * @param matrix Matrix to multiply with the scalar
         * @param scalar Scalar value to multiply
         * @return Resultant matrix after scalar multiplication
         */
        friend Matrix operator*(const Matrix& matrix, double scalar);

        /*!
         * @brief Scalar division operator
         *
         * @param matrix Matrix to divide
         * @param scalar Scalar value to divide by
         * @return Resultant matrix after scalar division
         */
        friend Matrix operator/(const Matrix& matrix, double scalar);

        /*!
         * @brief Scalar addition operator
         *
         * @param matrix Matrix to add
         * @param scalar Scalar value to add
         * @return Resultant matrix after scalar addition
         */
        friend Matrix operator+(const Matrix& matrix, double scalar);

        /*!
         * @brief Accessor operator for matrix elements
         *
         * @param i Row index
         * @param j Column index
         * @return Reference to the matrix element at (i, j)
         */
        double& operator()(int i, int j) const;

        //! Method to print the matrix
        void print();

        /*!
         * @brief Calculates the Euclidean norm of the matrix
         *
         * @return Euclidean norm of the matrix
         */
        double norm() const;

        /*!
         * @brief Calculates the cross product of the matrix with another matrix
         *
         * @param matrix2 Another matrix to compute the cross product with
         * @return Resultant matrix after cross product
         */
        Matrix cross(const Matrix& matrix2) const;

        //! Method to get the number of rows in the matrix
        int rows() const;

        //! Method to get the number of columns in the matrix
        int columns() const;

        /*!
         * @brief Calculates the transpose of the matrix
         *
         * @return Transpose of the matrix
         */
        Matrix traspuesta() const;

        /*!
         * @brief Calculates the inverse of the matrix
         *
         * @return Inverse of the matrix
         */
        Matrix inversa() const;

        /*!
         * @brief Extracts a submatrix from the matrix
         *
         * @param f Starting row index
         * @param l Ending row index
         * @param i Starting column index
         * @return Submatrix extracted from the matrix
         */
        Matrix subMatrix(int f, int l, int i) const;

        /*!
         * @brief Concatenates the matrix with another matrix
         *
         * @param matrix2 Another matrix to concatenate with
         * @return Concatenated matrix
         */
        Matrix concatenar(const Matrix& matrix2) const;

        /*!
         * @brief Checks if the matrix is equal to another matrix within a specified tolerance
         *
         * @param matrix2 Another matrix to compare with
         * @param TOL_ Tolerance value
         * @return True if matrices are equal within tolerance, false otherwise
         */
        bool isEqual(const Matrix& matrix2, double TOL_) const;

        /*!
         * @brief Generates an identity matrix of the specified size
         *
         * @param size Size of the identity matrix
         * @return Identity matrix of the specified size
         */
        static Matrix identidad(int size);


    private:
        //! Method to initialize the matrix.
        void initMatrix();

    private:
        //! Number of rows in the matrix
        int fil;
        //! Number of columns in the matrix
        int col;
        //! 2D array representing the matrix
        double **matrix;
};

#endif
