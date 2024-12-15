#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include "mathj.h"
#include "forceField.h"
#include <cstddef>

/// @brief Class for Molecule
/// @param _nAtoms Number of Atomes
/// @param _nBonds Number of Bonds
/// @param _nAngles Number of Angles
/// @param _nDihedrals Number of Dihedrals
/// @param _atoms Atom Number of Each Atoms
/// @param _bonds The Two Atoms that Form the Bonds
/// @param _angles The Three Atoms that Form the Angles
/// @param _dihedrals The Four Atoms that Form the Dihedrals
/// @param _pos Pointer to Matrix Object of Position, The Size Should Be nAtoms * DIMENSION
/// @param _bondLength Pointer to Matrix Object of Bond Lengths
/// @param _AngleDegree Pointer to Matrix Object of Angle Degrees
/// @param _DihedralDegree Pointer to Matrix Object Dihedral Degrees
/// @param _bondC Force Constant of Bonds
/// @param _angleC Force Constant of Angles
/// @param _dihedralC Force Constant of Dihedrals
class Molecule
{
private:
    size_t _nAtoms;
    size_t _nBonds;
    size_t _nAngles;
    size_t _nDihedrals;
    size_t _nPairs;

    unsigned int  *_atoms;
    unsigned int **_bonds;
    unsigned int **_angles;
    unsigned int **_dihedrals;
    unsigned int **_pairs;


    Matrix *_pos;
    Matrix *_interCoor;
    Matrix *_e;
    Matrix *_gq;
    Matrix *_gx;
    Matrix *_B;
    Matrix *_para1;
    Matrix *_para2;
    Matrix *_M;
    
    // Matrix *_bondLength;
    // Matrix *_angleDegree;
    // Matrix *_dihedralDegree;
    // Matrix *_bondC;
    // Matrix *_angleC;
    // Matrix *_dihedralC;
public:
    Molecule(const size_t& nAt, const size_t& nBo, unsigned int *at, unsigned int **bo, Matrix *pos);
    ~Molecule();
    void print();
    void output(std::ofstream &file, size_t format, double tol);

    /// @brief Generate Molecule Object Using Input File
    /// @param fileName The Input File
    /// @return Return Pointer to Molecule Object 
    static Molecule* readFromFile(const char* const fileName);
    void calInterCoord(const Matrix &posA, Matrix &interC, Matrix &B);
    void calInterCoord(const Matrix &posA, Matrix &interC);
    void findForceConst();
    void calEnergy(Matrix &e, Matrix &interC);
    void calGradient(Matrix &gq, Matrix &interC);
    void evolve(std::ofstream &file, unsigned int mode, double tol = 0.001);
    double lineSearch(std::ofstream &file, const Matrix& p, double a = 1, bool ifOutput = false);
};

#endif