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


    double calculateIC(const std::vector<Matrix> &pos);
    std::vector<double> calculateB(const std::vector<Matrix> &pos);

};


#endif