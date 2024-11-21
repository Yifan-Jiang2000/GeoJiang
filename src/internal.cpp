#include "internal.h"
#include <vector>
#include "const.h"
#include <cmath>
#include <iostream>
#include <iomanip>


Internal::Internal(Internal::InternalType t): _type(t) {}


double Internal::calculateIC(const std::vector<Matrix> &pos)
{
    double res;
    switch(_type)
    {
    case STRETCH:
    case VDW:
    {
        return norm(pos[0]-pos[1]);
    }
    case BEND:
    {
        Matrix rba = pos[0] - pos[1];
        Matrix rbc = pos[2] - pos[1];
        double theta = std::acos(dot(rba,rbc) / (norm(rba) * norm(rbc)));
        return theta;
    }
    case TORSION:
    {
        Matrix rab = pos[1] - pos[0];
        Matrix rbc = pos[2] - pos[1];
        Matrix rcd = pos[3] - pos[2];
        Matrix t = cross(rab, rbc);
        Matrix u = cross(rbc, rcd);
        Matrix v = cross(t, u);
        double cosphi = dot(t, u) / (norm(t) * norm(u));
        double sinphi = dot(rbc, v) / (norm(rbc) * norm(t) * norm(u));
        double phi = atan2(sinphi, cosphi);
        return phi;
    }
    default:
        std::cout << "Internal Coord Type NOT Defined!" << std::endl;
        exit(1);
    }
}

std::vector<double> Internal::calculateB(const std::vector<Matrix> &pos)
{
    std::vector<double> res;
    switch(_type)
    {
    case STRETCH:
    case VDW:
    {
        Matrix rba = pos[0] - pos[1];
        double rabL = norm(rba);
        Matrix da = (1.0 / rabL) * rba;
        for (size_t i = 0; i < DIMENSION; ++i) res.push_back(da._Mdata[i]);
        for (size_t i = 0; i < DIMENSION; ++i) res.push_back(-da._Mdata[i]);
        return res;
    }
    case BEND:
    {
        Matrix rba = pos[0] - pos[1];
        Matrix rbc = pos[2] - pos[1];
        double rabL = norm(rba);
        double rbcL = norm(rbc);
        Matrix p = cross(rba, rbc);
        double pL = norm(p);
        Matrix da = cross(rba, p) / (std::pow(rabL, 2) * pL);
        Matrix dc = cross(-rbc, p) / (std::pow(rbcL, 2) * pL);
        for (size_t i = 0; i < DIMENSION; ++i) res.push_back(da._Mdata[i]);
        for (size_t i = 0; i < DIMENSION; ++i) res.push_back(-da._Mdata[i] - dc._Mdata[i]);
        for (size_t i = 0; i < DIMENSION; ++i) res.push_back(dc._Mdata[i]);
        return res;
    }
    case TORSION:
    {
        Matrix rab = pos[1] - pos[0];
        Matrix rbc = pos[2] - pos[1];
        Matrix rcd = pos[3] - pos[2];
        Matrix rac = pos[2] - pos[0];
        Matrix rbd = pos[3] - pos[1];
        Matrix t = cross(rab, rbc);
        Matrix u = cross(rbc, rcd);
        double tL = norm(t);
        double uL = norm(u);
        double rbcL = norm(rbc);
        Matrix t1 = cross(t, rbc) / (std::pow(tL, 2) * rbcL);
        Matrix u1 = cross(u, rbc) / (std::pow(uL, 2) * rbcL);
        Matrix da = cross(t1, rbc);
        Matrix dd = cross(-u1, rbc);
        Matrix db = cross(rac, t1) + cross(-u1, rcd);
        Matrix dc = cross(t1, rab) + cross(rbd, -u1);
        for (size_t i = 0; i < DIMENSION; ++i) res.push_back(da._Mdata[i]);
        for (size_t i = 0; i < DIMENSION; ++i) res.push_back(db._Mdata[i]);
        for (size_t i = 0; i < DIMENSION; ++i) res.push_back(dc._Mdata[i]);
        for (size_t i = 0; i < DIMENSION; ++i) res.push_back(dd._Mdata[i]);
        return res;
    }
    default:
        std::cout << "Internal Coord Type NOT Defined!" << std::endl;
        exit(1);
    }
}