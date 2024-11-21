#ifndef HESSIAN_H
#define HESSIAN_H

#include "mathj.h"


class Hessian
{
public:
    static void guessM(Matrix &M);
    static void guessM(Matrix &M, size_t nB, size_t nA, size_t nD);
    static void updateM(Matrix &M, const Matrix &s, const Matrix &y);

};

#endif