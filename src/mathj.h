#ifndef _MATHJ_H
#define _MATHJ_H

#include <cstddef>
#include <fstream>
#include <cassert>


struct Matrix
{
    size_t _row;
    size_t _column;
    size_t _lD;
    double *_Mdata;
    bool _isSlice;

    Matrix();
    Matrix(size_t r, size_t c);
    Matrix(size_t r, size_t c, double *data);
    Matrix(const Matrix &m);
    Matrix(Matrix &&m);
    ~Matrix();
    
    void setZero();
    
    /// @brief Set to Identity Matrix, Only For Square Matrix
    void setIdentity();

    void setUnitVec();
    /// @brief Exchange Row a and Row b
    /// @param a a Starts with 0
    /// @param b b Starts with 0
    void swapRow(size_t a, size_t b);
    
    void swapColumn(size_t a, size_t b);
    
    //Subtract row a with k times row b

    /// @brief Subtract Row a with e Times of Row b
    /// @param a a Starts with 0
    /// @param b b Starts with 0
    /// @param e The Scaling Factor of Row b
    void subtractRow(size_t a, size_t b, double e);

    void divideRow(size_t a, double e);
 
    Matrix transpose() const;
    Matrix inverse() const;
    void House(Matrix &v, double &beta) const;


    /// @brief Find tridiagonal matrix out of a symmetric matrix using household transformation
    /// @param A Target matrix
    /// @param[out] V Orthogonal transformation matrix
    /// @param[out] W Working matrix
    /// @param[out] H Tridiagonal matrix
    void unblockedSymHessenberg(Matrix &A, Matrix &V, Matrix &W, Matrix &H) const;

    void unblockedHessenberg(Matrix &A, Matrix &T, bool ifSym = false) const;

    void QRDecompositionTri(Matrix &Q, Matrix &v);
    void QRDecomposition(Matrix &Q, Matrix &R) const;


    /// @brief Unblocked method for eigen decomposition
    /// @param[out] v Orthogonal transformation matrix (eigen vectors)
    void unblockedEigenDecomposition(Matrix &v);
    // void blockedEigenDecomposition(Matrix &v) const;
    void eigenDecomposition(Matrix &lambda, Matrix &v) const;

    /// @brief Call lapack dsyev routine to find eigenvalue and eigen vectors
    /// @param[out] lambda Eigen values
    /// @param[out] v Eigen vectors
    void eigenDecompositionLAPACK(Matrix &lambda, Matrix &v) const;
    Matrix getInvDiag() const;
    Matrix getInv() const;


    /// @brief Find one of the eigenvalues of a 2*2 matrix
    /// @return An eigenvalue
    double findEigen22() const;
    Matrix guassEliminaton() const;
    double det() const;
    double sum() const;
    void print() const;
    void print(std::ofstream &file) const;
    
    /**
    * @brief Cut few rows out of a Matrix, share the same memory. If want to cut columns, please take transpose first
    * @param start Start row, the first row is zero
    * @param end End row, the last row is _row - 1
    */
    Matrix cutR(size_t start, size_t end) const;
    Matrix cut(size_t startR, size_t endR, size_t startC, size_t endC) const;
    Matrix operator= (const Matrix &m);
    Matrix operator= (Matrix &&m);
    Matrix operator- ();
    void operator+= (const Matrix &m);
    void operator-= (const Matrix &m);
    void operator*= (const double &e);

    //Hadamard product
    void operator*= (const Matrix &m);
    void operator/= (const double &e);

    //Hadamard division
    void operator/= (const Matrix &m);

    inline double& operator() (size_t r, size_t c) const { 
        assert(r < _row && c < _column);
        return _Mdata[r * _lD + c]; 
    }
};

Matrix operator+ (const Matrix &l, const Matrix &r);

Matrix operator+ (const double &a, const Matrix &m);

Matrix operator- (const Matrix &l, const Matrix &r);

Matrix operator* (const Matrix &l, const Matrix &r);

Matrix operator* (const double &a, const Matrix &m);

Matrix operator/ (const Matrix &m, const double &s);


Matrix Hadamard(const Matrix &l, const Matrix &r);

Matrix blockMul(const Matrix &A, const Matrix &B, size_t blockSize, char blockType);
void selfBlockMul(Matrix &A, const Matrix &B, size_t blockSize, char blockType);

Matrix HadamardDiv(const Matrix &l, const Matrix &r);

Matrix pow(const Matrix &m, const double &s);

double dot(const Matrix &va, const Matrix &vb);
double norm(const Matrix &v);
Matrix cross(const Matrix &va, const Matrix &vb);
Matrix cos(const Matrix &m);
Matrix sin(const Matrix &m);


bool findInArray(unsigned int *ats, size_t size, unsigned int target);
inline int sign(double x) { return x < 0 ? -1 : 1; }
inline bool compareSize(const Matrix &a, const Matrix &b) { return a._row == b._row && a._column == b._column; }
// size_t permuteEigen(Matrix &D, Matrix &v, const Matrix &z);
void sortAdj(double *value, size_t *order, size_t size);

extern void (Matrix::*eigenSolver)(Matrix&, Matrix&) const;

#endif

