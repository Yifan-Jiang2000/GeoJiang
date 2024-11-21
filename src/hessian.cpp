#include "hessian.h"
#include <cassert>
#include <iostream>
#include <cmath>

void Hessian::guessM(Matrix &M)
{
    assert(M._row == M._column && "Hessian Matrix and its Inverse Should be Square Matrix!!");
    M.setZero();
    for (size_t i = 0; i < M._row; ++i) M(i, i) = 1.0 / 300.0;
    return;
}

void Hessian::guessM(Matrix &M, size_t nB, size_t nA, size_t nD)
{
    assert(M._row == M._column && "Hessian Matrix and its Inverse Should be Square Matrix!!");
    M.setZero();
    size_t i;
    for (i = 0; i < nB; ++i) M(i, i) = 1.0 / 600.0;
    for (;i < nB + nA; ++i) M(i, i) = 1.0 / 150.0;
    for (;i < nB + nA + nD; ++i) M(i, i) = 1.0 / 80.0;
}

void Hessian::updateM(Matrix &M, const Matrix &s, const Matrix &y)
{
    Matrix v = M * y;
    double dsy = dot(s, y), dyv = dot(y, v);
    M += ((dsy + dyv) / std::pow(dsy, 2)) * s * s.transpose() - (v * s.transpose() + s * v.transpose()) / dsy;
}


