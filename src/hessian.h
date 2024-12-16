#ifndef HESSIAN_H
#define HESSIAN_H

#include "mathj.h"


class Hessian
{
public:
    /// @brief Produce intial guess of inverse Hessian in cartesian coordiantes
    /// @param[out] M Inverse Hessian
    static void guessM(Matrix &M);

    /// @brief Produce intial guess of inverse Hessian in internal coordiantes
    /// @param M Inverse Hessian
    /// @param nB Number of bonds
    /// @param nA Number of angles
    /// @param nD Number of dihedrals

    static void guessM(Matrix &M, size_t nB, size_t nA, size_t nD);
    
    /// @brief Update inverse Hessian in internal coordiantes
    /// @param M Inverse Hessian
    /// @param s 
    /// @param y 
    static void updateM(Matrix &M, const Matrix &s, const Matrix &y);

};

#endif