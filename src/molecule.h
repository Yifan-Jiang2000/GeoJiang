#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include "mathj.h"
#include "forceField.h"
#include <cstddef>

/// @brief Class for Molecule
/// @param _nAtoms Number of atoms
/// @param _nBonds Number of bonds
/// @param _nAngles Number of angles
/// @param _nDihedrals Number of dihedrals
/// @param _nPairs Number of atom pairs

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

    /// @brief calculate internal coordinates of a molecule and the wilson B matrix
    /// @param posA cartesian coordiantes
    /// @param[out] interC internal coordinates
    /// @param[out] B wilson B matrix
    void calInterCoord(const Matrix &posA, Matrix &interC, Matrix &B);

    /// @brief calculate internal coordinates of a molecule
    /// @param posA cartesian coordiantes
    /// @param[out] interC internal coordinates 
    void calInterCoord(const Matrix &posA, Matrix &interC);

    /// @brief get a set of force parameters for each internal coordinates
    void findForceConst();

    /// @brief calculate energy from internal coordinates
    /// @param[out] e energy contribution of each internal coordinates
    /// @param interC internal coordinates
    void calEnergy(Matrix &e, Matrix &interC);

    /// @brief calculate energy gradient from internal coordinates
    /// @param[out] gq energy gradient contribution of each internal coordiantes
    /// @param interC internal coordinates
    void calGradient(Matrix &gq, Matrix &interC);


    /// @brief optimization process
    /// @param file output file
    /// @param mode optimization mode
    /// @param tol tolerance of gradient
    void evolve(std::ofstream &file, unsigned int mode, double tol = 0.001);

    /// @brief line search process for cartesian optimization
    /// @param file output file
    /// @param p predicted coordinate change
    /// @param a scale of the coordinate change
    /// @param ifOutput true in debugging mode, false in info only mode
    /// @return 
    double lineSearch(std::ofstream &file, const Matrix& p, double a = 1, bool ifOutput = false);
};

#endif