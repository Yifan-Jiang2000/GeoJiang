#ifndef ATOMTYPE_H
#define ATOMTYPE_H

#include <unordered_map>
#include <string>



/// @brief Convertion between atomic symbol and atomic number
class AtomType
{
private:
    static std::unordered_map<std::string, unsigned int> _symToNumber;
    static std::unordered_map<unsigned int, std::string> _numberToSym;
public:
    /// @brief Get Atom Number
    /// @param atomSym Atom Symbol
    /// @return Return the Atom Number
    static unsigned int lookupAtomNumber(const std::string& atomSym);

    /// @brief Get atomic symbol
    /// @param atomNum 
    /// @return Atomic symbol
    static std::string lookupAtomSym(const unsigned int& atomNum);
};

#endif