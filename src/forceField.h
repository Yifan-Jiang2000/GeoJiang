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
    std::array<double,2> findPara(std::vector<unsigned int> ats);

    void calEnergy(Matrix&& e, const Matrix& inter, const Matrix& para1, const Matrix& para2);
    void calGradient(Matrix&& g, const Matrix& inter, const Matrix& para1, const Matrix& para2);

};


#endif