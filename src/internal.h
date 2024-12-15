#ifndef INTERNAL_H
#define INTERNAL_H

#include <vector>
#include "mathj.h"

class Internal
{
public:
    enum InternalType
    {
        STRETCH,
        BEND,
        TORSION,
        VDW
    } _type;

    Internal(Internal::InternalType t);


    /// @brief Calculate internal coordinates using cartesian coordinates
    /// @param pos cartesian coordinates
    /// @return internal coordinates:bond length, angle or dihedral angle
    double calculateIC(const std::vector<Matrix> &pos);

    /// @brief Calculate element of wilson B matrix
    /// @param pos 
    /// @return Element of wilson B matrix
    std::vector<double> calculateB(const std::vector<Matrix> &pos);

};


#endif