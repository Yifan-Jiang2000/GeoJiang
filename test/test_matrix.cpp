#include "mathj.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <random>
#include <ctime>

const double TOL = 1e-8;

bool compareMatrix(const Matrix &c, const Matrix &d) {
    assert((c._row == d._row && c._column == d._column) && "The matrix should have the same size");
    for (size_t i = 0; i < c._row; ++i) {
        for (size_t j = 0; j < c._column; ++j)
            if (fabs(c(i, j) - d(i, j)) > TOL) {
                // std::cout << i << "   " << j << std::endl;
                return false;
            }
    }
    return true;
}



int main() {
    unsigned int res = 0;
    Matrix a = Matrix(20, 20);
    Matrix b = Matrix(20, 20);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dis_float(-1.0, 1.0);

    for (size_t i = 0; i < a._row * a._column; ++i)
    {
        a._Mdata[i] = dis_float(gen);
        b._Mdata[i] = dis_float(gen);
    }
    Matrix c = Matrix(20, 20);
    Matrix d = Matrix(20, 20);

    //Addition

    c = a + b;
    for (size_t i = 0; i < a._row * a._column; ++i) {
        d._Mdata[i] = a._Mdata[i] + b._Mdata[i];
    }
    if (compareMatrix(c, d)) {
        std::cout << "Addition of Matrix : Success" << std::endl;
    } else {
        ++res;
        std::cout << "Addition of Matrix : Fail" << std::endl;
    }

    //Subtraction
    c = a - b;
    for (size_t i = 0; i < a._row * a._column; ++i) {
        d._Mdata[i] = a._Mdata[i] - b._Mdata[i];
    }
    if (compareMatrix(c, d)) {
        std::cout << "Subtraction of Matrix : Success" << std::endl;
    } else {
        ++res;
        std::cout << "Subtraction of Matrix : Fail" << std::endl;
    }

    //Hadamard
    c = Hadamard(a, b);
    for (size_t i = 0; i < a._row * a._column; ++i) {
        d._Mdata[i] = a._Mdata[i] * b._Mdata[i];
    }
    if (compareMatrix(c, d)) {
        std::cout << "Hadamard Product of Matrix : Success" << std::endl;
    } else {
        ++res;
        std::cout << "Hadamard Product of Matrix : Fail" << std::endl;
    }    

    //Hadamard Division
    c = HadamardDiv(a, b);
    for (size_t i = 0; i < a._row * a._column; ++i) {
        d._Mdata[i] = a._Mdata[i] / b._Mdata[i];
    }
    if (compareMatrix(c, d)) {
        std::cout << "Hadamard Division of Matrix : Success" << std::endl;
    } else {
        ++res;
        std::cout << "Hadamard Division of Matrix : Fail" << std::endl;
    }    

    //Matrix Multiplication
    c = a * b;
    d.setZero();
    for (size_t i = 0; i < a._row; ++i) {
        for (size_t k = 0; k < b._row; ++k) {
            for (size_t j = 0; j < b._column; ++j) {
                d(i, j) += a(i, k) * b(k, j);
            }
        }
    }
    if (compareMatrix(c, d)) {
        std::cout << "Matrix Multiplication of Matrix : Success" << std::endl;
    } else {
        ++res;
        std::cout << "Matrix Multiplication of Matrix : Fail" << std::endl;
    }    


    //Transpose
    c = a.transpose();
    for (size_t i = 0; i < a._row; ++i) {
        for (size_t j = 0; j < a._column; ++j) {
            d(j, i) = a(i, j);
        }
    }
    if (compareMatrix(c, d)) {
        std::cout << "Transpose of Matrix : Success" << std::endl;
    } else {
        ++res;
        std::cout << "Transpose of Matrix : Fail" << std::endl;
    }    

    Matrix va(20, 1), vb(20, 1);
    for (size_t i = 0; i < va._row * va._column; ++i)
    {
        va._Mdata[i] = dis_float(gen);
        vb._Mdata[i] = dis_float(gen);
    }
    Matrix vtc(1, 20), vtd(1, 20);
    Matrix vc(20, 1), vd(20, 1);

    // Vector Transpose
    vtc = va.transpose();
    for (size_t i = 0; i < va._row; ++i) {
        for (size_t j = 0; j < va._column; ++j) {
            vtd(j, i) = va(i, j);
        }
    }

    if (compareMatrix(vtc, vtd)) {
        std::cout << "Transpose of Vector : Success" << std::endl;
    } else {
        ++res;
        std::cout << "Transpose of Vector : Fail" << std::endl;
    }        

    double dotc, dotd = 0;
    // Vector Dot Product
    dotc = dot(va, vb);
    dotd = (va.transpose() * vb)._Mdata[0];
    if (fabs(dotc - dotd) < TOL) {
        std::cout << "Dot Product of Vector : Success" << std::endl;
    } else {
        ++res;
        std::cout << "Dot Product of Vector : Fail" << std::endl;
    }

    // Matrix slice
    size_t startR = 5, endR = 14, startC = 5, endC = 14;
    Matrix sliceC = a.cut(startR, endR, startC, endC);
    Matrix sliceD(10, 10);
    for (size_t i = 0; i < sliceD._row; ++i) {
        for (size_t j = 0; j < sliceD._column; ++j) {
            sliceD(i, j) = a(i + 5, j + 5);
        }
    }
    if (compareMatrix(sliceC, sliceD)) {
        std::cout << "Slice of Matrix : Success" << std::endl;
    } else {
        ++res;
        std::cout << "Slice of Matrix : Fail" << std::endl;
    }        

    //Slice Addition
    Matrix SliceA = a.cut(startR, endR, startC, endC);
    Matrix SliceB = b.cut(startR, endR, startC, endC);
    Matrix C = SliceA + SliceB;
    Matrix DFull = a + b;
    Matrix D = DFull.cut(startR, endR, startC, endC);
    if (compareMatrix(C, D)) {
        std::cout << "Slice Addition of Matrix : Success" << std::endl;
    } else {
        ++res;
        std::cout << "Slice Addition of Matrix : Fail" << std::endl;
    }            
    // SliceA.print();
    // SliceB.print();

    //Block Multiplication
    Matrix blockI(20, 20);
    for (size_t i = 0; i < blockI._row * blockI._column; ++i)
    {
        blockI._Mdata[i] = dis_float(gen);
    }
    

    blockI.cut(0, 9, 0, 9).setIdentity();
    blockI.cut(10, 19, 0, 9).setZero();
    blockI.cut(0, 9, 10, 19).setZero();
    // blockI.print();
    c = a * blockI;
    d = blockMul(a, blockI, 10, 'I');
    if (compareMatrix(c, d)) {
        std::cout << "Block Multiplication of Matrix : Success" << std::endl;
    } else {
        ++res;
        std::cout << "Block Multiplication of Matrix : Fail" << std::endl;
    }                



    // Matrix v(a._row, 1);
    // double beta;
    // Matrix v_ = v.cut(2, a._row - 1, 1, 1);
    // a.cut(2, a._row - 1, 1, 1).House(v_, beta);

    //Hessenberg
    Matrix A = a, T = a;
    a.unblockedHessenberg(A, T);
    // A.print();
    bool isHessen = true;
    for (size_t i = 2; i < A._row; ++i) {
        for (size_t j = 0; j < i - 1; ++j) {
            if (fabs(A(i, j) > TOL)) isHessen = false;
        }
    }
    if (!isHessen) {
        std::cout << "Hessenberg of Matrix : Fail : Is not Hessenberg" << std::endl;
        ++res;
    } else {
        Matrix A_ = T.transpose() * a * T;
        if (compareMatrix(A, A_)) {
            std::cout << "Hessenberg of Matrix : Success" << std::endl;
        } else {
            ++res;
            std::cout << "Hessenberg of Matrix : Fail : Orthogonal Similarity Transformation is Wrong" << std::endl;
        }
    }

    // Symetric Hessenberg
    Matrix sym(10, 10);    
    for (size_t i = 0; i < sym._row; ++i) {
        for (size_t j = i; j < sym._column; ++j) {
            sym(i, j) = dis_float(gen);
            sym(j, i) = sym(i, j);
        }
    }

    Matrix V(sym._row, sym._column);
    Matrix W(sym._row, sym._column);
    Matrix Asym = sym;
    Matrix H(sym._row, sym._column);
    sym.unblockedSymHessenberg(Asym, V, W, H);
    // Asym.print();
    Matrix Asym_ = sym - V * W.transpose() - W * V.transpose();
    // Asym_.print();
    Matrix Asym__ = H.transpose() * sym * H;

    // Asym__.QRDecompositionTri(V, W);
    // Asym__.print();
    // (V.transpose() * Asym_ * V).print();


    // Unblocked Eigenvalue (Diagonalization) of Symetric Matrix
    Matrix lambda(sym._row, sym._column), v(sym._row, sym._column);
    v.setIdentity();
    Asym_.unblockedEigenDecomposition(v);
    // Asym_.print();
    bool isDiag = true;
    for (size_t i = 1; i < Asym_._row; ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (fabs(Asym_(i, j) > TOL)) isDiag = false;
        }
    }
    if (!isDiag) {
        std::cout << "Unblocked Tridiagnonal Eigenvalue of Matrix : Fail : Is not Diagonal" << std::endl;
        ++res;
    } else {
        Matrix lambda_ = v.transpose() * Asym__ * v;
        if (compareMatrix(Asym_, lambda_)) {
            std::cout << "Unblocked Tridiagnonal Eigenvalue of Matrix : Success" << std::endl;
        } else {
            ++res;
            std::cout << "Unblocked Tridiagnonal Eigenvalue of Matrix : Fail : Orthogonal Similarity Transformation is Wrong" << std::endl;
        }
    }


    // Eigenvalue of Symmetric Matrix
    v.setIdentity();
    sym.eigenDecomposition(lambda, v);
    isDiag = true;
    for (size_t i = 1; i < Asym_._row; ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (fabs(Asym_(i, j) > TOL)) isDiag = false;
        }
    }
    if (!isDiag) {
        std::cout << "Eigenvalue of Symmetric Matrix : Fail : Is not Diagonal" << std::endl;
        ++res;
    } else {
        Matrix lambda_ = v.transpose() * sym * v;
        lambda.print();
        if (compareMatrix(lambda, lambda_)) {
            std::cout << "Eigenvalue of Symmetric Matrix : Success" << std::endl;
        } else {
            ++res;
            std::cout << "Eigenvalue of Symmetric Matrix : Fail : Orthogonal Similarity Transformation is Wrong" << std::endl;
        }
    }

    //Sort Value and Order

    double value[20] = {-0.13,-0.81,0.36,-0.29,-0.99,0.76,0.04,0.7,0.84,-0.3,-0.21,0.43,-0.63,-0.53,-0.79,0.27,-0.85,0.73,0.21,-0.88};
    size_t order[20] = {   10,    9,  17,    3,   19,  11,  14,  7,  16,   8,   13,   2,    6,    5,    0,  12,   15,   4,   1,   18};
    bool success = true;
    sortAdj(value, order, 20);
    // std::cout << "!!!" << std::endl;
    // for (size_t i = 0; i < 20; ++i) {
    //     std::cout << value[i] << " " << order[i] << std::endl;
    // }
    for (size_t i = 0; i < 19; ++i) {
        if (value[i] > value[i + 1]) {
            success = false;
            break;
        }
    }        
    if (success) {
        std::cout << "Sort Value and Order : Success" << std::endl;
    } else {
        ++res;
        std::cout << "Sort Value and Order : Fail" << std::endl;
    }


 
    // Permutation of Eigenvalue

    // lambda.swapColumn(1, 2);
    // lambda.swapRow(1, 2);
    // v.swapColumn(1, 2);
    // Matrix lambda_ = v.transpose() * sym * v;
    // // lambda.print();
    // // lambda_.print();
    // if (compareMatrix(lambda, lambda_)) {
    //     std::cout << "Permutation of Eigenvalue : Success" << std::endl;
    // } else {
    //     ++res;
    //     std::cout << "Permutation of Eigenvalue : Fail" << std::endl;
    // }




    // sym.print();
    // lambda.print();
    // (v.transpose() * sym * v).print();
    // Asym__.print();

    // Matrix Q(sym._row, sym._column);
    // Asym__.QRDecompositionTri(Q);
    // Q.print();
    // Asym__.print();
    // (Asym__ * Q).print();

    // Matrix ASym = sym;
    // Matrix TSym = sym;
    // sym.unblockedHessenberg(ASym, TSym, true);
    // ASym.print();
    // for (size_t i = 2; i < ASym._row; ++i) {
    //     for (size_t j = 0; j < i - 1; ++j) {
    //         if (fabs(ASym(i, j) > TOL)) isHessen = false;
    //     }
    // }
    // if (!isHessen) {
    //     std::cout << "Symmetric Hessenberg of Matrix : Fail : Is not Hessenberg" << std::endl;
    //     ++res;
    // } else {
    //     Matrix A_ = TSym.transpose() * sym * TSym;
    //     if (compareMatrix(ASym, A_)) {
    //         std::cout << "Symmetric Hessenberg of Matrix : Success" << std::endl;
    //     } else {
    //         ++res;
    //         std::cout << "Symmetric Hessenberg of Matrix : Fail : Orthogonal Similarity Transformation is Wrong" << std::endl;
    //     }       
    // }


    // 500 * 500 eigenvalue time
    // A.print();
    // Matrix symBIG(200, 200);
    // Matrix vBIG(symBIG._row, symBIG._column);
    // Matrix lambdaBIG(symBIG._row, symBIG._column);
    // // Matrix WBIG(symBIG._row, symBIG._column);
    // // Matrix HBIG(symBIG._row, symBIG._column);
    // for (size_t i = 0; i < symBIG._row; ++i) {
    //     for (size_t j = i; j < symBIG._column; ++j) {
    //         symBIG(i, j) = dis_float(gen);
    //         symBIG(j, i) = symBIG(i, j);
    //     }
    // }
    // std::clock_t start = clock();
    // // symBIG.unblockedSymHessenberg()
    // symBIG.eigenDecomposition(lambdaBIG, vBIG);
    // // symBIG.eigenDecomposition(lambdaBIG, vBIG);
    // std::clock_t end = clock();
    // std::cout << "200 * 200 Eigen Decomposition : "<< (double)(end - start) / CLOCKS_PER_SEC << " s" << std::endl;
    // start = clock();
    // // symBIG.unblockedSymHessenberg()
    // symBIG.eigenDecompositionLAPACK(lambdaBIG, vBIG);
    // // symBIG.eigenDecomposition(lambdaBIG, vBIG);
    // end = clock();
    // std::cout << "200 * 200 Eigen Decomposition(Lapack) : "<< (double)(end - start) / CLOCKS_PER_SEC << " s" << std::endl;// Matrix AsBIG = symBIG, TsBIG = symBIG;
    // symBIG.unblockedSymHessenberg(AsBIG, VBIG, WBIG, TsBIG);
    // for (size_t i = 0; i < symBIG._row; ++i) {
    //     for (size_t j = i + 2; j < symBIG._column; ++j) {
    //         if (fabs(AsBIG(i, j)) > TOL || fabs(AsBIG(j, i)) > TOL) {
    //             std::cout << i << "  " << j << "  " << AsBIG(j, i) << std::endl;
    //         }
    //     }
    // }
    // As.print();

    // A.print();
    // A_.print();
    
    // T.setIdentity();
    // A = a;
    // Matrix A__ = A.cut(12, 19, 12, 19);
    // Matrix T__ = T.cut(12, 19, 12, 19);
    // a.cut(12,19,12,19).blockHessenberg(A__, T__);
    // // A.print();
    // A = T.transpose() * a * T;
    // // A.print();

    // Matrix A___ = A.cut(4, 11, 4, 11);
    // Matrix T___ = T.cut(4, 11, 4, 11);
    // A.cut(4, 11, 4, 11).blockHessenberg(A___, T___);
    // A.print();
    // A = T.transpose() * a * T;
    // A.print();

    return res;
}