#ifndef PARAMETER_H
#define PARAMETER_H

#include <vector>
#include "const.h"
#include "mathj.h"
#include <array>


template<size_t n>
struct ForceConst
{
    unsigned int _atoms[n];
    double _p1;
    double _p2;


    /// @brief compare atoms with those corresponding to the force constant
    /// @param ats 
    /// @return if atoms match
    bool ifAtomMatch(std::vector<unsigned int> ats);
};



class MolForcePara
{
public:
    enum ForceType
    {
        fSTRETCH,
        fBEND,
        fTORSION,
        fVDW
    } _type;
    


    static std::vector<ForceConst<N_ATOM_STRETCH>> strecthPara;
    static std::vector<ForceConst<N_ATOM_BEND>> bendPara;
    static std::vector<ForceConst<N_ATOM_TORSION>> torsionPara;
    static std::vector<ForceConst<N_ATOM_VDW>> vdwPara;


    MolForcePara(MolForcePara::ForceType t);


    static int ifSet;
    static void setPara();
    /// @brief find the corresponding parameter
    /// @param ats 
    /// @return the two parameter
    std::array<double,2> findPara(std::vector<unsigned int> ats);

    /// @brief calculate energy
    /// @param[out] e energy
    /// @param inter internal coordinates
    /// @param para1 the first parameter
    /// @param para2 the second parameter
    void calEnergy(Matrix&& e, const Matrix& inter, const Matrix& para1, const Matrix& para2);
    /// @brief calculate gradient
    /// @param[out] g gradient
    /// @param inter internal coordiantes
    /// @param para1 the first parameter
    /// @param para2 the second parameter
    void calGradient(Matrix&& g, const Matrix& inter, const Matrix& para1, const Matrix& para2);

};


#endif