#include "mathj.h"
#include <iostream>
#include <cstring>
#include <cmath>
#include <cstddef>
#include <cassert>
#include <iomanip>
#include "f2c.h"
#include "clapack.h"

#define EPSILON 1e-12


void (Matrix::*eigenSolver)(Matrix&, Matrix&) const = &Matrix::eigenDecompositionLAPACK;


Matrix::Matrix():_row(0),_column(0), _Mdata(nullptr), _isSlice(false){}

Matrix::Matrix(size_t r, size_t c):_row(r),_column(c), _lD(c), _isSlice(false)
{
    _Mdata = (double*)malloc(r * c * sizeof(double));
}

Matrix::Matrix(size_t r, size_t c,double *data):_row(r),_column(c), _lD(c), _isSlice(false)
{
    _Mdata = (double*)malloc(r * c * sizeof(double));
    for (size_t i = 0; i < r * c; ++i) {
        _Mdata[i] = data[i];
    }
}

Matrix::Matrix(const Matrix &m):_row(m._row),_column(m._column), _lD(m._column), _isSlice(false)
{
    _Mdata = (double*)malloc(_row * _column * sizeof(double));
    for (size_t i = 0; i < _row; ++i) {
        for (size_t j = 0; j < _column; ++j) {
            (*this)(i, j) = m(i, j);
        }
    }
}

Matrix::Matrix(Matrix &&m):_row(m._row),_column(m._column), _lD(m._lD), _Mdata(m._Mdata), _isSlice(m._isSlice)
{
    m._Mdata = nullptr;
}

void Matrix::setZero()
{
    for (size_t i = 0; i < _row; ++i) {
        for (size_t j = 0; j < _column; ++j) {
            (*this)(i, j) = 0;
        }
    }
    return;
}

void Matrix::setIdentity()
{
    if (_row != _column)
        std::cout << "It is not a square matrix!!" << std::endl;
    setZero();
    for (size_t i = 0; i < _row; ++i)
        (*this)(i, i) = 1;
    return;
}

void Matrix::setUnitVec()
{
    assert((_row == 1 || _column == 1) && "Vector should only have 1 row or 1 column.");
    for (size_t i = 0; i < _row; ++i) {
        for (size_t j = 0; j < _column; ++j) {
            (*this)(i, j) = 1;
        }
    }
    return;
}

Matrix::~Matrix()
{
    if (!_isSlice)
        free(_Mdata);
}

void Matrix::swapRow(size_t a, size_t b)
{
    double temp;
    for (size_t i = 0; i < _column; ++i) {
        temp = (*this)(a, i);
        (*this)(a, i) = (*this)(b, i);
        (*this)(b, i) = temp;
    }
    // double *temp = (double*)malloc(_column * sizeof(double));
    // memcpy(temp, _Mdata + a * _lD, _column * sizeof(double));
    // memcpy(_Mdata + a * _lD, _Mdata + b * _column, _column * sizeof(double));
    // memcpy(_Mdata + b * _lD, temp, _column * sizeof(double));
    // free(temp);
}

void Matrix::swapColumn(size_t a, size_t b) {
    double temp;
    for (size_t i = 0; i < _row; ++i) {
        temp = (*this)(i, a);
        (*this)(i, a) = (*this)(i, b);
        (*this)(i, b) = temp;
    }
}

void Matrix::subtractRow(size_t a, size_t b, double e)
{
    for (int j = 0; j < _column; ++j)
        (*this)(a, j) -= e * (*this)(b, j);
        // _Mdata[a * _column + j] -= e * _Mdata[b * _column + j];
}

void Matrix::divideRow(size_t a, double e)
{
    for (int j = 0; j < _column; ++j)
        (*this)(a, j) /= e;
        // _Mdata[a * _column + j] /= e;
    return;
}

Matrix Matrix::transpose() const
{
    Matrix tMat = Matrix(_column, _row);
    // if (_row == 1 || _column == 1) {
    //     memcpy(tMat._Mdata, _Mdata, _row * _column * sizeof(double));
    // } else {
    for (size_t i = 0; i < _row; ++i) {
        for (size_t j = 0; j < _column; ++j) {
            tMat(j, i) = (*this)(i, j);
        }
    }
    // tMat._Mdata[j * _row + i] = _Mdata[i * _column + j];
    // }
    return tMat;
}

Matrix Matrix::inverse() const
{
    assert(_row == _column && "It is not a square matrix!!");
    // if (_row != _column)
    //     std::cout << "It is not a square matrix!!" << std::endl;
    Matrix iMat= Matrix(_column, _row);
    Matrix temp = *this;
    iMat.setIdentity();
    for (size_t j = 0; j < temp._column - 1; ++j)
    {   
        if (fabs(temp._Mdata[j * temp._column + j]) < EPSILON)
        {
            for (size_t p = j + 1; p < temp._row; ++p)
            {
                if (fabs(temp._Mdata[p * temp._column + j]) > EPSILON)
                {
                    temp.swapRow(j, p);
                    iMat.swapRow(j, p);
                    break;
                }
            }
        }
        if (fabs(temp._Mdata[j * temp._column + j]) < EPSILON)
            continue;
        for (size_t i = j + 1; i < temp._row; ++i)
        {
            double e = temp(i, j) / temp(j, j);
            temp.subtractRow(i, j, e);
            iMat.subtractRow(i, j, e);
        }
    }
    // temp.print();
    for (size_t i = 0; i < temp._row; ++i)
    {
        double diag = temp._Mdata[i * temp._column + i];
        if (fabs(diag) < EPSILON)
            std::cout << "This matrix do not have a inverse" << std::endl;
        for (size_t j = 0; j < temp._column; ++j)
        {
            temp(i, j) /= diag;
            iMat(i, j) /= diag;
        }
    }
    for (size_t j = temp._column - 1; j > 0; --j)
    {
        for (size_t i = 0; i < j; ++i)
        {
            double e = temp(i, j);
            temp.subtractRow(i, j, e);
            iMat.subtractRow(i, j, e);
        }
    }
    // temp.print();
    return iMat;
}

void Matrix::House(Matrix &v, double &beta) const {
//     double mu = dot(*this, *this);
//     double sigma = mu - (*this)(0, 0) - (*this)(0, 0);
//     v = *this;
//     v(0, 0) = 1.;
//     if (sigma < EPSILON) {
//         beta = 0.;
//         return;
//     }
//     mu = std::pow(mu, 0.5);
//     if ((*this)(0, 0) > 0) {
//         v(0, 0) = -sigma / ((*this)(0, 0) + mu);
//     } else {
//         v(0, 0) = (*this)(0, 0) - mu;
//     }
//     beta = 2 * std::pow(v(0, 0), 2) / (sigma + std::pow(v(0, 0), 2));
//     v /= v(0, 0);
    // std::cout << beta << std::endl;
    // v.print();
    // return;
    v = *this;
    v(0, 0) += (double)sign(v(0, 0)) * norm(v);
    beta = 2 / dot(v, v);
    return;
}

void Matrix::unblockedSymHessenberg(Matrix &A, Matrix &V, Matrix &W, Matrix &H) const {
    V.setZero();
    W.setZero();
    A = *this;
    H.setIdentity();
    // Matrix Adummy = A;
    for (size_t i = 0; i < _column - 2; ++i) {
        // if (i > 0) {
        //     Matrix A__ = Adummy.cut(i, A._row - 1, i, i);
        //     Matrix vvT = V.cut(i, i, 0, i - 1).transpose();
        //     Matrix wwT = W.cut(i, i, 0, i - 1).transpose();
        //     Matrix vv = V.cut(i, V._row - 1, 0, i - 1);
        //     Matrix ww = W.cut(i, W._row - 1, 0, i - 1);
        //     A__ -= vv * wwT + ww * vvT;
        // }
        // A.print();
        // Adummy.print();
        Matrix v_ = V.cut(i + 1, V._row - 1, i, i);
        Matrix h_ = H.cut(0, H._row - 1, i + 1, H._column - 1);
        double beta;
        A.cut(i + 1, A._row - 1, i, i).House(v_, beta);
        // double beta_;
        // Matrix v__ = v_;
        // Adummy.cut(i + 1, A._row - 1, i, i).House(v__, beta_);
        // v_.print();
        // v__.print();
        // std::cout << beta << "  " << beta_ << std::endl;
        h_ -= beta * h_ * v_ * v_.transpose();
        v_._Mdata -= v_._lD;
        v_._row += 1;
        Matrix w_ = W.cut(i, W._row - 1, i, i);
        Matrix p = beta * A.cut(i, _row - 1, i, _row - 1) * v_;
        w_ = p - (beta * dot(v_, p) / 2) * v_;
        Matrix A_ = A.cut(i, A._row - 1, i, A._column - 1);
        A_ -= v_ * w_.transpose() + w_ * v_.transpose();
    }
    return;
}



void Matrix::unblockedHessenberg(Matrix &A, Matrix &T, bool ifSym) const {
    assert(compareSize(*this, A));
    assert(compareSize(*this, T));
    A = *this;
    T.setIdentity();
    Matrix H(_row, _column);
    for (size_t i = 0; i < _column - 2; ++i) {
        H.setIdentity();
        // A.print();
        // Matrix a = A.cut(i + 1, _row - 1, i, i);
        // a.print();
        // if (fabs(a(0, 0)) < EPSILON) continue;
        double beta;
        Matrix v(_row, 1);
        v.setZero();
        Matrix v_ = v.cut(i + 1, _row - 1, 0, 0);
        
        A.cut(i + 1, _row - 1, i, i).House(v_, beta);
        // Matrix b(_row - i - 1, 1);
        // b.setZero();
        // b(0, 0) = 1;
        // Matrix u = a + (double)sign(a(0, 0)) * norm(a) * b;
        // u /= norm(u);
        Matrix Hr = H.cut(i + 1, _row - 1, i + 1, _column - 1);
        // v_.print();
        // v_.transpose().print();
        // Hr.print();
        // std::cout << v._lD << std::endl;
        v_ * v_.transpose();
        // std::cout << "!!" << std::endl;
        Hr -= beta * v_ * v_.transpose();
        Matrix HT = H.transpose();
        // Matrix w = beta * A * v;
        if (ifSym) {
            Matrix A_ = A.cut(i, _row - 1, i, _column - 1);
            Matrix v__ = v.cut(i, _row - 1, 0, 0);
            Matrix p = beta * (A_ * v__);
            Matrix w = p - (beta * dot(v__, p) / 2) * v__;
            A_ -= (w * v__.transpose() + v__ * w.transpose());
            A.print();
            if (i == 2) {
                // v.print();
                // w.print();
                // A.print();
            }
        } else {
            A = H * A * HT;
            // A.print();
        }
        // if (i == 2) T.print();        
        // if (i == 2) {
            // Matrix T_ = T;
            // Matrix T__ = T_.cut(0, T_._row - 1, i + 1, T_._column - 1);
            // T__ -= beta * T__ * v_ * v_.transpose();
            // T_.print();
        // }
        selfBlockMul(T, HT, i + 1, 'I');
        // if (i == 2) T.print();

        // T = T * HT;
    }
    return;
}

void Matrix::blockHessenberg(Matrix &A, Matrix &T) const {
    if (_row < 32) {
        unblockedHessenberg(A, T);
        return;
    } else {
        A = *this;
        T.setIdentity();
        Matrix A_ = A, T_ = T;
    
    }
}

void Matrix::getUpHessenberg(Matrix &trans, Matrix &hessen) const
{
    hessen = *this;
    trans.setIdentity();
    Matrix H(_row, _column);
    for (size_t i = 0; i < hessen._column - 2; ++i)
    {
        // std::clock_t start = clock();
        Matrix a(hessen._row - i - 1, 1);
        // memcpy(a._Mdata, hessen._Mdata + i * hessen._column + i + 1, (hessen._column - i - 1) * sizeof(double));
        for (size_t j = i + 1; j < hessen._row; ++j)
        {
            a._Mdata[j - i - 1] = hessen._Mdata[j * hessen._column + i];
        }
        // std::cout << a._Mdata[0] << std::endl;
        Matrix b(hessen._row - i - 1, 1);
        b.setZero();
        b._Mdata[0] = 1;        
        Matrix u = a + (double)sign(a._Mdata[0]) * norm(a) * b;
        if (fabs(a._Mdata[0]) < EPSILON) continue;
        Matrix n = u / norm(u);
        Matrix I(hessen._row - i - 1, hessen._column - i - 1);
        I.setIdentity();
        Matrix Hr = I - 2 * n * n.transpose();
        H.setZero();
        for (size_t m = 0; m < Hr._row; ++m)
        {
            // memcpy(H._Mdata + (m + i + 1) * H._column + i + 1, Hr._Mdata + m * Hr._column, Hr._column - 1 * sizeof(double));
            for (size_t n = 0; n < Hr._column; ++n)
            {
                H._Mdata[(m + i + 1) * H._column + n + i + 1] = Hr._Mdata[m * Hr._column + n];
            }
        }
        for (size_t m = 0; m < H._row - Hr._row; ++m)
        {
            H._Mdata[m *  H._column + m] = 1;
        }
        // H.print();
        // std::cout << "!!!" << std::endl;
        hessen = H * hessen * H.transpose();
        trans = trans * H.transpose();
        // std::cout << "!!!" << std::endl;
        for (size_t m = 0; m < hessen._row * hessen._column; ++m)
            hessen._Mdata[m] = fabs(hessen._Mdata[m]) < EPSILON ? 0 : hessen._Mdata[m];
        // std::clock_t end = clock();
        // std::cout << (double)(end - start) / CLOCKS_PER_SEC << std::endl;
    }
    // (trans.transpose() * hessen * trans).print();
    // hessen.print();
    // this->print();
    return;
}


void Matrix::QRDecompositionTri(Matrix &Q, Matrix &v)
{
    Q.setIdentity();
    Matrix u(2, 1);
    double beta;
    // Matrix R = *this;
    for (size_t i = 0; i < _column - 1; ++i) {
        this->cut(i, i + 1, i, i).House(u, beta);
        Matrix uuT = beta * u * u.transpose();
        Matrix R_ = this->cut(i, i + 1, i, _column - 1);
        R_ -= uuT * R_;
        Matrix Q_ = Q.cut(0, Q._row - 1, i, i + 1);
        Q_ -= Q_ * uuT;
        Matrix v_ = v.cut(0, v._row - 1, i, i + 1);
        v_ -= v_ * uuT;
        // Matrix thisCut = this->cut(i, i + 1, i, _column - 1);
        // thisCut -= uuT * thisCut;
        // Matrix thisCut_ = this->cut(0, Q._row - 1, i, i + 1);
        // thisCut_ -= thisCut_ * uuT;
    }
    // R.print();
    // *this = R * Q;
    return;
}


void Matrix::QRDecomposition(Matrix &Q, Matrix &R) const
{
    assert(_row == _column && "Matrix should be square to do QR Decomposition");
    Q.setIdentity();
    Matrix P_new(_row, _column);
    R = *this;
    for (size_t i = 0; i < R._row - 1; ++i)
    {
        // R.print();
        Matrix a(R._column - i, 1);
        for (size_t j = i; j < R._row; ++j)
        {
            a._Mdata[j - i] = R._Mdata[j * _column + i];
        }
        // a.print();
        if (fabs(a._Mdata[0]) < EPSILON) continue;
        Matrix b(R._column - i, 1);
        b.setZero();
        b._Mdata[0] = 1;
        Matrix u = a + (double)sign(a._Mdata[0]) * norm(a) * b;
        Matrix n = u / norm(u);
        Matrix I(R._column - i, R._column - i);
        I.setIdentity();
        Matrix P = I - 2 * n * n.transpose();
        P_new.setZero();
        for (size_t m = 0; m < P._row; ++m)
        {
            for (size_t n = 0; n < P._column; ++n)
            {
                P_new._Mdata[(m + i) * P_new._column + n + i] = P._Mdata[m * P._column + n];
            }
        }
        for (size_t m = 0; m < P_new._row - P._row; ++m)
        {
            P_new._Mdata[m *  P_new._column + m] = 1;
        }
        R = P_new * R;
        Q = Q * P_new;
        for (size_t m = 0; m < R._row * R._column; ++m)
            R._Mdata[m] = fabs(R._Mdata[m]) < EPSILON ? 0 : R._Mdata[m];
        // R.print();
    }
    // (R * Q).print();
}

void Matrix::unblockedEigenDecomposition(Matrix &v) {
    // v.setIdentity();
    const double tol = 1e-9;
    Matrix Q(this->_row, this->_column);
    // Matrix R(this->_row, this->_column);
    // Matrix I(this->_row, this->_column);
    Matrix Q_ = Q.cut(0, Q._row - 1, 0, Q._column - 1);
    Matrix R_ = this->cut(0, this->_row - 1, 0, this->_column - 1);
    // Matrix I_ = I.cut(0, I._row - 1, 0, I._column - 1);
    Matrix v_ = v.cut(0, v._row - 1, 0, v._column - 1);
    // I.setIdentity();        
    // bool isConverge;
    // unsigned int notZero = 0;
    unsigned int cycles = 0;
    double sqrt2 = std::pow(2, 0.5);
    double mu;
    do
    {
        // if (R_(2, 1) < sqrt2 * R_(1, 0)) {
        Matrix m22 = R_.cut(0, 1, 0, 1);
        mu = m22.findEigen22();
        // } else {
            // mu = R_(0, 0);
        // }
        for (size_t i = 0; i < R_._row; ++i) R_(i, i) -= mu;
        R_.QRDecompositionTri(Q_, v_);
        R_ = R_ * Q_;
        for (size_t i = 0; i < R_._row; ++i) R_(i, i) += mu;
        // isConverge = true;
        // if (fabs(R_(1, 0)) < tol && R_._row == 2) break;
        while (R_._row > 1 && fabs(R_(1, 0)) < tol) {
            // std::cout << R_(1, 0) << std::endl;
            // std::cout << R_._row << std::endl;
            // this->print();
            --R_._row; --R_._column; R_._Mdata += R_._lD + 1;
            --Q_._row; --Q_._column; Q_._Mdata += Q_._lD + 1;
            --v_._column; v_._Mdata += 1;
            // if (R_._row == 2) cycles = 0;
        }
        cycles++;
        if (cycles > 1000000) {
            this->print();
            break;
        }
        // std::cout << R_._row << std::endl;
    } while(R_._row > 1);
    // std::cout << cycles << std::endl;
    size_t *order = new size_t[_row];
    double *diag = new double[_row];
    for (size_t i = 0; i < _row; ++i) {
        order[i] = i;
        diag[i] = (*this)(i, i);
    }
    
    sortAdj(diag, order, _row);

    Matrix vtemp(_row, _column);
    for (size_t i = 0; i < _column; ++i) {
        for (size_t j = 0; j < _row; ++j) {
            vtemp(j, i) = v(j, order[i]);
        }
    }

    for (size_t i = 0; i < _row; ++i) {
        (*this)(i, i) = diag[i];
    }

    v = vtemp;
    return;
    // this->print();
}


#if 0
void Matrix::blockedEigenDecomposition(Matrix &v)
{
    v.setIdentity();
    const unsigned int maxBlockSize = 25;
    // Matrix W(this->_row, this->_column), U(this->_row, this->_column);
    // this->unblockedSymHessenberg(lambda, U, W, v);
    if (this->row > maxBlockSize) {
        unsigned int bisect = int(_row / 2)
        double alpha = (*this)(blockSize, blockSize - 1);
        (*this)(blockSize - 1, blockSize - 1) -= alpha;
        (*this)(blockSize - 1, blockSize) -= alpha;
        (*this)(blockSize, blockSize - 1) -= alpha;
        (*this)(blockSize, blockSize) -= alpha;

        Matrix lambda1 = this->cut(0, blockSize - 1, 0, blockSize - 1);
        Matrix lambda2 = this->cut(blockSize, _row - 1, blockSize, _column - 1);
        Matrix v1 = v.cut(0, blockSize - 1, 0, blockSize - 1);
        Matrix v2 = v.cut(blockSize, v._row - 1, blockSize, v._column - 1);

        lambda1.blockedEigenDecomposition(v1);
        lambda2.blockedEigenDecomposition(v2);

        Matrix z(_row, 1);
        Matrix z1 = z.cut(0, blockSize - 1, )


    } else {
        lambda.unblockedEigenDecomposition(v);
        // std::cout << cycles << std::endl;
    }

}
#endif

void Matrix::eigenDecomposition(Matrix &lambda, Matrix &v) const {
    v.setIdentity();
    Matrix W(_row, _column), V(_row, _column);
    this->unblockedSymHessenberg(lambda, V, W, v);
    lambda.unblockedEigenDecomposition(v);
}

void Matrix::eigenDecompositionLAPACK(Matrix &lambda, Matrix &v) const
{
    v = *this;
    double *v_ = new double[v._row * v._column];
    memcpy(v_, v._Mdata, v._row * v._column * sizeof(double));
    long N = v._row, LDA = v._row;
    double *eigen = new double[v._row];
    long LWORK = 4 * v._row;
    double *WORK = new double[LWORK];
    long INFO;
    dsyev_("V", "U", &N, v_, &LDA, eigen, WORK, &LWORK, &INFO);
    lambda.setZero();
    for (size_t i = 0; i < v._row; ++i) lambda._Mdata[i * lambda._column + i] = eigen[i];
    memcpy(v._Mdata, v_, v._row * v._column * sizeof(double));
    v = v.transpose();
}


Matrix Matrix::getInvDiag() const
{
    assert(_row == _column && "getInvDiag:The matrix should be square matrix");
    Matrix temp = *this;
    for (size_t i = 0; i < temp._row; ++i)
    {
        temp._Mdata[i * _column + i] = fabs(temp._Mdata[i * _column + i]) > EPSILON ? 1 / temp._Mdata[i * _column + i] : 0.0;
    }
    return temp;
}

Matrix Matrix::getInv() const
{
    assert(_row == _column && "getInv:The matrix should be square matrix");
    Matrix lambda(_row, _column);
    Matrix v(_row, _column);
    (this->*eigenSolver)(lambda, v);
    Matrix invLambda = lambda.getInvDiag();
    return v * invLambda * v.transpose();
}

double Matrix::findEigen22() const
{
    assert((_row == 2 && _column == 2) && "findEigen22 : Matrix should be 2 * 2!!");
    double m = ((*this)(0, 0) + (*this)(1, 1)) / 2;
    double p = ((*this)(0, 0) * (*this)(1, 1)) - ((*this)(0, 1) * (*this)(1, 0));
    double msq = pow(m, 2);
    if (msq < p)
        return m;
    double deviation = pow(msq - p, 0.5);
    double a = m + deviation;
    double b = m - deviation;
    return fabs(a) - (*this)(1, 1) > fabs(b) - (*this)(1, 1) ? b : a;
}

double Matrix::powerIter(Matrix &v)
{
    const double tol = 1e-5;
    Matrix vec(_row, 1);
    vec.setUnitVec();
    Matrix vec_old = vec;
    double sigma = 0, sigma_old;
    do
    {
        sigma_old = sigma;
        vec = (*this) * vec_old;
        sigma = vec._Mdata[0] / vec_old._Mdata[0];
        vec_old = vec / norm(vec);
        // std::cout << vec._Mdata[0] << std::endl;
        std::cout << sigma << std::endl;
    } while (fabs(sigma - sigma_old) > tol);
    Matrix temp = *this;
    for (size_t i = 0; i < temp._row; ++i) temp._Mdata[i * temp._column + i] -= sigma;
    temp = temp.guassEliminaton();
    for (size_t i = 0; i < temp._row; ++i) v._Mdata[i] = temp._Mdata[i * _column + _column - 1];
    v._Mdata[v._row - 1] = 1;

    return sigma;
}

Matrix Matrix::guassEliminaton() const
{
    assert(_row == _column);
    Matrix temp = *this;
    unsigned int swapCount = 0;
    for (size_t i = 0; i < temp._row; ++i)
    {
        if (fabs(temp._Mdata[i * temp._column + i]) < EPSILON)
        {
            for (size_t j = i + 1; j < temp._row; ++j)
            {
                if (fabs(temp._Mdata[j * temp._column + i]) > EPSILON)
                {
                    temp.swapRow(i, j);
                    ++swapCount;
                }
            }
        }
        if (fabs(temp._Mdata[i * temp._column + i]) < EPSILON) continue;
        temp.divideRow(i, temp._Mdata[i * temp._column + i]);
        for (size_t j = 0; j < temp._row; ++j)
        {
            if (i == j) continue;
            temp.subtractRow(j, i, temp._Mdata[j * temp._column + i]);
        }
    }
    if (swapCount % 2 == 1)
        temp = -temp;
    return temp;
}

double Matrix::det() const
{
    if (_row != _column)
        std::cout << "It is not a square matrix!!" << std::endl;
    Matrix temp = *this;
    unsigned int swapCount = 0;
    for (size_t j = 0; j < temp._column - 1; ++j)
    {   
        if (fabs(temp._Mdata[j * temp._column + j]) < EPSILON)
        {
            for (size_t p = j + 1; p < temp._row; ++p)
            {
                if (fabs(temp._Mdata[p * temp._column + j]) > EPSILON)
                {
                    temp.swapRow(j, p);
                    ++swapCount;
                    break;
                }
            }
        }
        if (fabs(temp._Mdata[j * temp._column + j]) < EPSILON)
            continue;
        for (size_t i = j + 1; i < temp._row; ++i)
        {
            double e = temp._Mdata[i * temp._column + j] / temp._Mdata[j * temp._column + j];
            temp.subtractRow(i, j, e);
        }
    }
    double res = 1;
    for (size_t i = 0; i < temp._row; ++i)
        res *= temp._Mdata[i * temp._column + i];
    if (swapCount % 2 != 0)
        res *= -1;
    return res;
}

double Matrix::sum() const
{
    double res = 0;
    for (size_t i = 0; i < _row * _column; ++i) res += _Mdata[i];
    return res;
}

void Matrix::print() const
{
    for (size_t i = 0; i < _row; ++i)
    {
        for (size_t j = 0; j < _column; ++j)
        {
            std::cout << std::setprecision(6) << std::fixed << std::setw(15) << (*this)(i, j) << "    ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Matrix::print(std::ofstream &file) const
{
    for (size_t i = 0; i < _row; ++i)
    {
        for (size_t j = 0; j < _column; ++j)
        {
            file << std::setprecision(8) << std::fixed << std::setw(15) << (*this)(i, j) << "    ";
        }
        file << std::endl;
    }
    file << std::endl;
}

Matrix Matrix::cutR(size_t start, size_t end) const {
    Matrix m(end - start + 1, _column);
    m._lD = _column;
    m._Mdata = _Mdata + start * _column;
    m._isSlice = true;
    return m;
}

Matrix Matrix::cut(size_t startR, size_t endR, size_t startC, size_t endC) const {
    Matrix m(endR - startR + 1, endC - startC + 1);
    m._lD = _lD;
    m._Mdata = _Mdata + startR * _lD + startC;
    m._isSlice = true;
    return m;
}

Matrix Matrix::operator= (const Matrix &m)
{
    // if (this == &m) return *this;
    assert((_row == m._row && _column == m._column) && "The Matrix should have exactly same size!!");
    for (size_t i = 0; i < m._row; ++i) {
        for (size_t j = 0; j < m._column; ++j) {
            (*this)(i, j) = m(i, j);
        }
    }
    return *this;
}

Matrix Matrix::operator= (Matrix &&m)
{
    if (this == &m) return *this;
    assert((_row == m._row && _column == m._column) && "The Matrix should have exactly same size!!");
    if (_isSlice) {
        for (size_t i = 0; i < m._row; ++i)
            for (size_t j = 0; j < m._column; ++j)
                (*this)(i, j) = m(i, j);
        return *this;
    } else {
        free(_Mdata);
        _Mdata = m._Mdata;
        m._Mdata = nullptr;
        return *this;
    }
}


Matrix Matrix::operator-()
{
    Matrix res(_row, _column);
    for (size_t i = 0; i < _row; ++i) {
        for (size_t j = 0; j < _column; ++j) {
            res(i, j) = - (*this)(i, j);
        }
    }
    return res;
}

void Matrix::operator+= (const Matrix &m)
{
    assert((_row == m._row && _column == m._column) && "The Matrix should have exactly same size!!");
    for (size_t i = 0; i < _row; ++i) {
        for (size_t j = 0; j < _column; ++j) {
            (*this)(i, j) += m(i, j);
        }
    }
    return;
}

void Matrix::operator-= (const Matrix &m)
{
    assert((_row == m._row && _column == m._column) && "The Matrix should have exactly same size!!");
    
    for (size_t i = 0; i < _row; ++i) {
        for (size_t j = 0; j < _column; ++j) {
            (*this)(i, j) -= m(i, j);
        }
    }
    return;
}

void Matrix::operator*= (const double &e)
{
    for (size_t i = 0; i < _row; ++i) {
        for (size_t j = 0; j < _column; ++j) {
            (*this)(i, j) *= e;
        }
    }
    return;
}

void Matrix::operator*= (const Matrix &m)
{
    assert((_row == m._row && _column == m._column) && "The Matrix should have exactly same size!!");
    for (size_t i = 0; i < _row; ++i) {
        for (size_t j = 0; j < _column; ++j) {
            (*this)(i, j) *= m(i, j);
        }
    }
}

void Matrix::operator/= (const double &e)
{
    for (size_t i = 0; i < _row; ++i) {
        for (size_t j = 0; j < _column; ++j) {
            (*this)(i, j) /= e;
        }
    }
    return;
}


void Matrix::operator/= (const Matrix &m)
{
    assert((_row == m._row && _column == m._column) && "The Matrix should have exactly same size!!");
    for (size_t i = 0; i < _row; ++i) {
        for (size_t j = 0; j < _column; ++j) {
            (*this)(i, j) /= m(i, j);
        }
    }
    return;
}

Matrix operator* (const Matrix &l, const Matrix &r)
{
    if (l._column != r._row) {
        std::cout << l._column << " != " << r._row << std::endl;
    }
    assert(l._column == r._row && "The column of the left matrix does not equal to the row of the right matrix!!");
    Matrix res = Matrix(l._row,r._column);
    res.setZero();
    // double aaaaaaa;
    for (size_t i = 0; i < l._row; ++i) {
        // std::cout << "i = " << i << std::endl;
        for (size_t k = 0; k < r._row; ++k) {
            // aaaaaaa = 0;
            // std::cout << "k = " << k << std::endl;
            for (size_t j = 0; j < r._column; ++j) {
                // std::cout << "j = " << j << std::endl;   
                res(i, j) += l(i, k) * r(k, j);
            }
        }
    }
    return res;
}

Matrix operator* (const double &a, const Matrix &m)
{
    Matrix res = m;
    res *= a;
    return res;
}

Matrix operator+ (const Matrix &l, const Matrix &r)
{
    assert((l._row == r._row && l._column == r._column) && "The Matrix should have exactly same size!!");
    Matrix res = l;
    res += r;
    return res;
}

Matrix operator+ (const double &a, const Matrix &m)
{
    Matrix res = Matrix(m._row,m._column);
    for (size_t i = 0; i < m._row; ++i) {
        for (size_t j = 0; j < m._column; ++j) {
            res(i, j) = m(i, j) + a;
        }
    }
    return res;
}

Matrix operator- (const Matrix &l, const Matrix &r)
{
    assert((l._row == r._row && l._column == r._column) && "The Matrix should have exactly same size!!");
    Matrix res = l;
    res -= r;
    return res;
}

Matrix operator/ (const Matrix &m, const double &s)
{
    assert(s != 0 && "The Denominator cannot be 0");
    Matrix res = m;
    res /= s;
    return res;
}

Matrix blockMul(const Matrix &A, const Matrix &B, size_t blockSize, char blockType) {
    Matrix res = A;
    if (blockType == 'I') {
        Matrix B22 = B.cut(blockSize, res._row - 1, blockSize, res._column - 1);
        Matrix res12 = res.cut(0, blockSize - 1, blockSize, res._column - 1);
        Matrix res22 = res.cut(blockSize, res._row - 1, blockSize, res._column - 1);
        res12 = res12 * B22;
        res22 = res22 * B22;
    }
    return res;
}

void selfBlockMul(Matrix &A, const Matrix &B, size_t blockSize, char blockType) {
    if (blockType == 'I') {
        Matrix B22 = B.cut(blockSize, B._row - 1, blockSize, B._column - 1);
        Matrix A12 = A.cut(0, blockSize - 1, blockSize, A._column - 1);
        Matrix A22 = A.cut(blockSize, A._row - 1, blockSize, A._column - 1);
        A12 = A12 * B22;
        A22 = A22 * B22;
    }
    return;
}




Matrix Hadamard(const Matrix &l, const Matrix &r)
{
    assert((l._row == r._row && l._column == r._column) && "The Matrix should have exactly same size!!");
    Matrix res = l;
    res *= r;
    return res;
}

Matrix HadamardDiv(const Matrix &l, const Matrix &r)
{
    assert((l._row == r._row && l._column == r._column) && "The Matrix should have exactly same size!!");
    Matrix res = l;
    res /= r;
    return res;
}

Matrix pow(const Matrix &m, const double &s)
{
    Matrix res = Matrix(m._row, m._column);
    for (size_t i = 0; i < m._row * m._column; ++i){
        for (size_t j = 0; j < m._column; ++j) {
            res(i, j) = std::pow(m(i, j), s);
        }
    }
    return res;
}

double dot(const Matrix &va, const Matrix &vb)
{
    assert((va._column == 1 && vb._column == 1) && "The Matrix in dot product should be of 1 column(vector)");
    assert(va._row == vb._row && "The Matrix in dot product should have same number of rows");
    double res = 0.0;
    for (size_t i = 0; i < va._row; ++i) res += va(i, 0) * vb(i, 0);
    return res;
}

double norm(const Matrix &v)
{
    assert(v._column == 1 && "The Matrix should be of 1 column(vector)");
    double res = 0.0;
    for (size_t i = 0; i < v._row; ++i) res += pow(v(i, 0), 2);
    res = pow(res, 0.5);
    return res;
}

Matrix cross(const Matrix &va, const Matrix &vb)
{
    assert((va._column == 1 && vb._column == 1) && "The Matrix in cross product should be of 1 column(vector)");
    assert((va._row == 3 && vb._row == 3) && "The Matrix in cross product should be of 3 rows");
    Matrix res(3,1);
    res(0, 0) = va(1, 0) * vb(2, 0) - va(2, 0) * vb(1, 0);
    res(1, 0) = va(2, 0) * vb(0, 0) - va(0, 0) * vb(2, 0);
    res(2, 0) = va(0, 0) * vb(1, 0) - va(1, 0) * vb(0, 0);
    return res;
}

Matrix cos(const Matrix &m)
{
    Matrix res(m._row, m._column);
    for (size_t i = 0; i < m._row; ++i) {
        for (size_t j = 0; j < m._column; ++j) {
            res(i, j) = cos(m(i, j));
        }
    }
    return res;
}

Matrix sin(const Matrix &m)
{
    Matrix res(m._row, m._column);    
    for (size_t i = 0; i < m._row; ++i) {
        for (size_t j = 0; j < m._column; ++j) {
            res(i, j) = sin(m(i, j));
        }
    }
    return res;
}


bool findInArray(unsigned int *ats, size_t size, unsigned int target)
{
    for (size_t i = 0; i < size; ++i) {
        if (ats[i] == target) return true;
    }
    return false;
}

// size_t permuteEigen(Matrix &D, Matrix &v, Matrix &z) {
//     size_t *order = new size_t[D._row];
//     double *diag = new double[D._row];
//     size_t sizeFlaten = 0;
//     for (size_t i = 0; i < D._row; ++i) {
//         order[i] = i;
//     }
//     for (size_t i = 0; i < D._row; ++i) {
//         diag[i] = D(i, i);
//     }
//     for (size_t i = 0; i < z._row; ++i) {
//         if (fabs(z(i, 0)) < EPSILON && i != sizeFlaten) {
//             double temp = z(i, 0);
//             z(i, 0) = z(i, 0);
//             z(i, 0) = temp;
//             temp = diag[i];
//             diag[i] = diag[sizeFlaten];
//             diag[sizeFlaten] = temp;
//             size_t otemp = order[i];
//             order[i] = order[sizeFlaten];
//             order[sizeFlaten] = otemp;
//             ++sizeFlaten;
//         }
//     }
//     sortAdj(diag + sizeFlaten, order + sizeFlaten, D._row - sizeFlaten);
//     Matrix v_(v._row, v._column);
//     for (size_t i = 0; i < v._column; ++i) {
//         Matrix v__ = v_.cut(0, v_._row - 1, i, i);
//         v__ = v._cut(0, v_._row, order[i], order[i]);
//     }
//     // v = v_;
//     Matrix z_(z._row, z._column);
//     for (size_t i = 0; i < z._row ++i) {
//         z_(i, 0) = z(order[i], 0);
//     }

//     for (size_t i = sizeFlaten; i < D._row - 1; ++i) {
//         size_t j = i;
//         while ( j < D._row && fabs(diag[j] - diag[j + 1] < EPSILON)) {
//             ++j;
//         }
//         --j;
//         if (j > i) {
//             Matrix z__ = z_.cutR(i, j);
//             double beta;
//             Matrix u(j - i + 1, 1);
//             z__.House(u, beta);
//             Matrix v__ = v_.cut(0, v_._row - 1, i, j);
//             v__ -= beta * v__ * u * u.transpose();
//             Matrix vtemp = v__;
//             for (size_t m = i; m > sizeFlaten - 1; --m) {
                
//             }
//         }
//     }
// }

void sortAdj(double *value, size_t *order, size_t size) {
    if (size < 2) return;
    size_t i = 0, j = size - 1;
    double pivot = value[0];
    size_t pivotOrder = order[0];
    bool right = true;
    if (value[j] < pivot) {
        value[i] = value[j];
        order[i] = order[j];
        right = false;
    }
    do {
        if(right) {
            --j;
            if (value[j] < pivot) {
                value[i] = value[j];
                order[i] = order[j];
                right = false;
            }
        } else {
            ++i;
            if (value[i] > pivot) {
                value[j] = value[i];
                order[j] = order[i];
                right = true;
            }
        }
    } while(i != j);
    value[i] = pivot;
    order[i] = pivotOrder;
    
    // std::cout << "!!!" << std::endl;
    if (i > 1){
        sortAdj(value, order, i);
    }
    if (i < size - 1){
        sortAdj(value + i + 1, order + i + 1, size - i - 1);
    }
    return;
}

