#include "atomType.h"
#include <unordered_map>
#include <iostream>

std::unordered_map<std::string, unsigned int> AtomType::_symToNumber = {
    {"H", 1},
    {"C", 6}
};

std::unordered_map<unsigned int, std::string> AtomType::_numberToSym = {
    {1, "H"},
    {6, "C"}
};

unsigned int AtomType::lookupAtomNumber(const std::string& atomSym)
{
    if (_symToNumber.find(atomSym) != _symToNumber.end()) 
        return _symToNumber[atomSym];
    else
    {
        std::cout << "The atom " << atomSym << " is NOT identified!" << std::endl;
        exit(1);
    }
}

std::string AtomType::lookupAtomSym(const unsigned int& atomNum)
{
    if (_numberToSym.find(atomNum) != _numberToSym.end()) 
        return _numberToSym[atomNum];
    else
    {
        std::cout << "The atom " << atomNum << " is NOT identified!" << std::endl;
        exit(1);
    }
}