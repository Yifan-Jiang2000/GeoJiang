#include <vector>
#include "forceField.h"
#include <unordered_map>
#include <cmath>
#include <array>
#include <iostream>


std::vector<ForceConst<N_ATOM_STRETCH>> MolForcePara::strecthPara;
std::vector<ForceConst<N_ATOM_BEND>> MolForcePara::bendPara;
std::vector<ForceConst<N_ATOM_TORSION>> MolForcePara::torsionPara;
std::vector<ForceConst<N_ATOM_VDW>> MolForcePara::vdwPara;
int MolForcePara::ifSet = 0;



template<size_t n>
bool ForceConst<n>::ifAtomMatch(std::vector<unsigned int> ats)
{
    bool isEqual = true;
    for (size_t i = 0; i < n; ++i)
    {
        if (ats[i] != _atoms[i])
        {
            isEqual = false;
            break;
        }
    }
    if (isEqual)
        return isEqual;
    isEqual = true;
    for (size_t i = 0; i < n; ++i)
    {
        if (ats[i] != _atoms[n - i - 1])
        {
            isEqual = false;
            break;
        }
    }
    return isEqual;
}


MolForcePara::MolForcePara(MolForcePara::ForceType t):_type(t) { if (ifSet == 0) setPara(); }

void MolForcePara::setPara()
{
    strecthPara.push_back({{1,6}, 350.0, 1.11});
    strecthPara.push_back({{6,6}, 300.0, 1.53});

    bendPara.push_back({{1,6,1}, 35.0, 109.5 * RADION});
    bendPara.push_back({{1,6,6}, 35.0, 109.5 * RADION});
    bendPara.push_back({{6,6,6}, 60.0, 109.5 * RADION});


    torsionPara.push_back({{1,6,6,1}, 0.300, 3});
    torsionPara.push_back({{1,6,6,6}, 0.300, 3});
    torsionPara.push_back({{6,6,6,6}, 0.300, 3});

    std::unordered_map<unsigned int, std::pair<double,double>> vdwAtom = {
        {1, {0.03, 1.2}},
        {6, {0.07, 1.75}}
    };

    for (const std::pair<unsigned int, std::pair<double,double>> &a:vdwAtom)
        for (const std::pair<unsigned int, std::pair<double,double>> &b:vdwAtom)
            if (a.first <= b.first)
                vdwPara.push_back({{a.first, b.first}, pow(a.second.first * b.second.first, 0.5), 2 * pow(a.second.second * b.second.second, 0.5)});

    ifSet = 1;
}

std::array<double,2> MolForcePara::findPara(std::vector<unsigned int> ats)
{
    switch(_type)
    {
    case fSTRETCH:
        for (auto c:strecthPara)
        {
            if (c.ifAtomMatch(ats))
                return std::array<double,2>({c._p1, c._p2});
        }
        std::cout << "Force Const NOT Found" << std::endl;
        exit(1);
    case fBEND:
        for (auto c:bendPara)
        {
            if (c.ifAtomMatch(ats))
                return std::array<double,2>({c._p1, c._p2});
        }
        std::cout << "Force Const NOT Found" << std::endl;
        exit(1);
    case fTORSION:
        for (auto c:torsionPara)
        {
            if (c.ifAtomMatch(ats))
                return std::array<double,2>({c._p1, c._p2});
        }
        std::cout << "Force Const NOT Found" << std::endl;
        exit(1);
    case fVDW:
        for (auto c:vdwPara)
        {
            if (c.ifAtomMatch(ats))
                return std::array<double,2>({c._p1, c._p2});
        }
        // std::cout << ats[0] << ats[1] << std::endl;
        std::cout << "Force Const NOT Found" << std::endl;
        exit(1);
    }
}


void MolForcePara::calEnergy(Matrix&& e, const Matrix& inter, const Matrix& para1, const Matrix& para2)
{
    switch(_type)
    {
    case fSTRETCH:
    case fBEND:
        e = Hadamard(para1, pow(inter - para2, 2));
        break;
    case fTORSION:
        e = Hadamard(para1,(1.0 + cos(Hadamard(para2, inter))));
        break;
    case fVDW:
        Matrix div = HadamardDiv(para2, inter);
        e = 4 * Hadamard(para1, (pow(div, 12)-pow(div, 6)));
        break;
    }

}

void MolForcePara::calGradient(Matrix&& g, const Matrix& inter, const Matrix& para1, const Matrix& para2)
{
    switch(_type)
    {
    case fSTRETCH:
    case fBEND:
        g = 2 * Hadamard(para1, inter - para2);
        break;
    case fTORSION:
        g = -Hadamard(para2, Hadamard(para1, sin(Hadamard(para2, inter))));
        break;
    case fVDW:
        g = 4 * 6 * Hadamard(para1, HadamardDiv(pow(para2, 6), pow(inter, 7)) - HadamardDiv(2 * pow(para2, 12), pow(inter, 13)));
        break;
    }
}
