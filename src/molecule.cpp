#include "molecule.h"
#include "const.h"
#include <string>
#include <iostream>
#include <cassert>
#include "atomType.h"
#include "mathj.h"
#include <vector>
#include <array>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "internal.h"
#include "hessian.h"
#include "f2c.h"
#include "clapack.h"
#include "blaswrap.h"
#include <cstring>





Molecule::Molecule(const size_t& nAt, const size_t& nBo, unsigned int *at, unsigned int **bo, Matrix *pos)
    :_nAtoms(nAt), _nBonds(nBo), _atoms(at), _bonds(bo), _pos(pos)
{

    _nPairs = _nAtoms * (_nAtoms - 1) / 2;
    _pairs = new unsigned int*[_nPairs];
    for (size_t i = 0; i < _nPairs; ++i)
    {
        _pairs[i] = new unsigned int[N_ATOM_VDW];
    }
    size_t count = 0;
    for (size_t i = 0; i < _nAtoms; ++i)
    {
        for (size_t j = i + 1; j < _nAtoms; ++j)
        {
            _pairs[count][0] = i + 1; 
            _pairs[count][1] = j + 1;
            ++count;
        }
    }

    // Get the Angles
    size_t nAn = 0;
    std::vector<std::array<size_t, N_ATOM_BEND>> angles;
    for (size_t i = 0; i < _nAtoms; ++i)
    {
        if (_atoms[i] == 1)
            continue;
        std::vector<size_t> bondedAtoms;
        for (size_t j = 0; j < _nBonds; ++j)
        {
            if (_bonds[j][0] == i + 1)
                bondedAtoms.push_back(_bonds[j][1]);

            if (_bonds[j][1] == i + 1)
                bondedAtoms.push_back(_bonds[j][0]);
        }
        for (size_t j = 0; j < bondedAtoms.size(); ++j)
        {
            for (size_t k = j + 1; k < bondedAtoms.size(); ++k)
            {
                ++nAn;
                angles.push_back({bondedAtoms[j], i + 1, bondedAtoms[k]});
            }
        }
    }
    _nAngles = nAn;
    _angles = new unsigned int*[_nAngles];
    for (size_t i = 0; i < _nAngles; ++i)
    {
        _angles[i] = new unsigned int[N_ATOM_BEND];
        for (size_t j = 0; j < N_ATOM_BEND; ++j)
            _angles[i][j] = angles[i][j];
    }

    //Get the Dihedrals

    size_t nDi = 0;
    std::vector<std::array<size_t, N_ATOM_TORSION>> dihedrals;

    for (size_t i = 0; i < _nBonds; ++i)
    {
        if (_atoms[_bonds[i][0] - 1] == 1 || _atoms[_bonds[i][1] - 1] == 1)
            continue;
        std::vector <size_t> bondedAtomA, bondedAtomB;
        for (size_t j = 0; j < _nBonds; ++j)
        {
            if (i == j)
                continue;
            if (_bonds[j][0] == _bonds[i][0])
                bondedAtomA.push_back(_bonds[j][1]);
            else if (_bonds[j][0] == _bonds[i][1])
                bondedAtomB.push_back(_bonds[j][1]);
            else if (_bonds[j][1] == _bonds[i][0])
                bondedAtomA.push_back(_bonds[j][0]);
            else if (_bonds[j][1] == _bonds[i][1])
                bondedAtomB.push_back(_bonds[j][0]);            
        }
        // if (_bonds[i][0] == 4 && _bonds[i][1] == 6)
        // {
        //     for (const auto& a:bondedAtomA)
        //     {
        //         std::cout << a << " ";
        //     }
        //     std::cout << std::endl;
        //     for (const auto& a:bondedAtomB)
        //     {
        //         std::cout << a << " ";
        //     }
        //     std::cout << std::endl;
        // }
        for (const auto& a:bondedAtomA)
        {
            for (const auto& b:bondedAtomB)
            {
                ++nDi;
                dihedrals.push_back({a, _bonds[i][0], _bonds[i][1], b});
            }
        }
    }
    _nDihedrals = nDi;
    _dihedrals = new unsigned int*[_nDihedrals];
    for (size_t i = 0; i < _nDihedrals; ++i)
    {
        _dihedrals[i] = new unsigned int[N_ATOM_TORSION];
        for (size_t j = 0; j < N_ATOM_TORSION; ++j)
            _dihedrals[i][j] = dihedrals[i][j];
    }

    //Create Matrix for Internal Coordinates
    _interCoor = new Matrix(_nBonds + _nAngles + _nDihedrals + _nPairs, 1);
    _gq = new Matrix(_interCoor->_row, 1);
    _e = new Matrix(_interCoor->_row, 1);
    _gx = new Matrix(DIMENSION * _nAtoms, 1);
    _B = new Matrix(_gq->_row, _gx->_row);
    calInterCoord(*_pos, *_interCoor, *_B);
    findForceConst();
    calEnergy(*_e, *_interCoor);
    calGradient(*_gq, *_interCoor);
    *(_gx) = _B->transpose() * *(_gq);

    _M = nullptr;
}

Molecule::~Molecule()
{
    delete[] _atoms;

    for (size_t i = 0; i < _nBonds; ++i)
        delete[] _bonds[i];
    delete[] _bonds;

    for (size_t i = 0; i < _nAngles; ++i)
        delete[] _angles[i];
    delete[] _angles;

    for (size_t i =0; i < _nDihedrals; ++i)
        delete[] _dihedrals[i];
    delete[] _dihedrals;

    for (size_t i =0; i < _nPairs; ++i)
        delete[] _pairs[i];
    delete[] _pairs;
    
    delete _pos;
    delete _interCoor;
    delete _para1;
    delete _para2;
    delete _e;
    delete _gq;
    delete _gx;
    delete _B;
    delete _M;
}

void Molecule::print()
{
    std::cout << "Positions" << std::endl;
    for (size_t i = 0; i < _nAtoms; ++i)
    {
        std::cout << _atoms[i] << " " << _pos->_Mdata[i] << " " << _pos->_Mdata[i + _nAtoms] << " " << _pos->_Mdata[i + 2 * +_nAtoms] << std::endl;
    }
    std::cout << "Bonds = " << _nBonds << std::endl;   
    for (size_t i = 0; i < _nBonds; ++i)
    {
        std::cout << _bonds[i][0] << " " << _bonds[i][1] << std::endl;
    }
    std::cout << "Angles = " << _nAngles << std::endl;   
    for (size_t i = 0; i < _nAngles; ++i)
    {
        std::cout << _angles[i][0] << " " << _angles[i][1] << " " << _angles[i][2] << std::endl;
    }
    std::cout << "Dihedrals = " << _nDihedrals << std::endl;   
    for (size_t i = 0; i < _nDihedrals; ++i)
    {
        std::cout << _dihedrals[i][0] << " " << _dihedrals[i][1] << " " << _dihedrals[i][2] << " " << _dihedrals[i][3] << std::endl;
    }
}

void Molecule::output(std::ofstream &file, size_t format, double tol)
{
    file << "-------------------------------------" << std::endl;
    file << "|This program is made by Yifan Jiang|" << std::endl;
    file << "|    This output is in format " << format << "     |" << std::endl;
    file << "-------------------------------------" << std::endl << std::endl;
        file << "The input file has: " << _nAtoms << " atoms" << std::endl;
        file << "Atoms and coordinates (in Å):" << std::endl;
        for (size_t i = 0; i < _nAtoms; ++i)
            file << AtomType::lookupAtomSym(_atoms[i]) << std::setprecision(6) << std::fixed <<
            std::setw(14) << _pos->_Mdata[i] <<
            std::setw(14) << _pos->_Mdata[i + _nAtoms] <<
            std::setw(14) << _pos->_Mdata[i + _nAtoms * 2] << std::endl;
        file << "Number of coordinates:" << std::endl;
        file << "Stretching:" << std::setw(7) << _nBonds <<
                " Bending:" << std::setw(7) << _nAngles <<
                " Torsion:" << std::setw(7) << _nDihedrals << std::endl;
        file << "Internal:" << std::setw(7) << _nBonds + _nAngles + _nDihedrals <<
                " Cartesian:" << std::setw(7) << DIMENSION * _nAtoms << " " << std::endl;
        file << "Potential energy at input structure:" << std::endl;
        file << std::setw(12) << std::setprecision(6) << std::fixed <<
                _e->sum() << " kcal/mol" << std::endl;
        file << "Stretch, Bend, Torsion, VDW components of potential energy:" << std::endl;
        file << std::setprecision(6) << std::fixed <<
                std::setw(12) << _e->cutR(0, _nBonds-1).sum() <<
                std::setw(13) << _e->cutR(_nBonds, _nBonds + _nAngles - 1).sum() <<
                std::setw(13) << _e->cutR(_nBonds + _nAngles, _nBonds + _nAngles + _nDihedrals - 1).sum() <<
                std::setw(13) << _e->cutR(_nBonds + _nAngles + _nDihedrals, _nBonds + _nAngles + _nDihedrals + _nPairs - 1).sum() <<
                std::endl;
    switch(format)
    {
    /*
        1 : singlepoint energy
        2 : cartesian optimization with debugging output
        3 : internal optimization with debugging output
        4 : cartesian optimization with less output
        6 : internal optimization with less output
    */
    case 1:
    {

        file << "List of all bonds: (At1 - At2, with labels, and distance in Angstrom, energy contrib in kcal/mol)" << std::endl;
        for (size_t i = 0; i < _nBonds; ++i)
            file << AtomType::lookupAtomSym(_atoms[_bonds[i][0] - 1]) << std::setw(4) << _bonds[i][0] <<
                    std::setw(4) << "-" <<
                    std::setw(4) << AtomType::lookupAtomSym(_atoms[_bonds[i][1] - 1]) << std::setw(4) << _bonds[i][1] << ":" <<
                    std::setprecision(5) << std::fixed <<
                    std::setw(13) << _interCoor->_Mdata[i] << 
                    std::setw(14) << _e->_Mdata[i] << std::endl;
        file << "List of all bending angles: (At1 - At2 - At3, with labels, angle in radian then degrees, energy contribution)" << std::endl;
        for (size_t i = 0; i < _nAngles; ++i)
            file << AtomType::lookupAtomSym(_atoms[_angles[i][0] - 1]) << std::setw(4) << _angles[i][0] <<
                    std::setw(4) << "-" <<
                    std::setw(4) << AtomType::lookupAtomSym(_atoms[_angles[i][1] - 1]) << std::setw(4) << _angles[i][1] <<
                    std::setw(4) << "-" << 
                    std::setw(4) << AtomType::lookupAtomSym(_atoms[_angles[i][2] - 1]) << std::setw(4) << _angles[i][2] << ":" <<
                    std::setprecision(6) << std::fixed <<
                    std::setw(13) << _interCoor->_Mdata[i + _nBonds] <<
                    std::setprecision(3) << 
                    std::setw(11) << _interCoor->_Mdata[i + _nBonds] * 180 / PI <<
                    std::setprecision(5) << 
                    std::setw(15) << _e->_Mdata[i + _nBonds] << std::endl;
        file << "List of all torsional angles: (At1 - At2 - At3 - At4, with labels, angle in radian then degrees,  energy contrib in kcal/mol)" << std::endl;
        for (size_t i = 0; i < _nDihedrals; ++i)
            file << AtomType::lookupAtomSym(_atoms[_dihedrals[i][0] - 1]) << std::setw(4) << _dihedrals[i][0] <<
                    std::setw(4) << "-" <<
                    std::setw(4) << AtomType::lookupAtomSym(_atoms[_dihedrals[i][1] - 1]) << std::setw(4) << _dihedrals[i][1] <<
                    std::setw(4) << "-" <<                    
                    std::setw(4) << AtomType::lookupAtomSym(_atoms[_dihedrals[i][2] - 1]) << std::setw(4) << _dihedrals[i][2] <<
                    std::setw(4) << "-" <<  
                    std::setw(4) << AtomType::lookupAtomSym(_atoms[_dihedrals[i][3] - 1]) << std::setw(4) << _dihedrals[i][3] << ":" <<
                    std::setprecision(6) << std::fixed <<
                    std::setw(13) << _interCoor->_Mdata[i + _nBonds + _nAngles] <<
                    std::setprecision(3) << 
                    std::setw(11) << _interCoor->_Mdata[i + _nBonds + _nAngles] * 180 / PI <<
                    std::setprecision(6) << 
                    std::setw(14) << _e->_Mdata[i + _nBonds + _nAngles] << std::endl;        
        file << "List of all unique atom pairs: (At1 - At2, with labels, distance,  vdW energy contrib in kcal/mol)" << std::endl;
        for (size_t i = 0; i < _nPairs; ++i)
            file << AtomType::lookupAtomSym(_atoms[_pairs[i][0] - 1]) << std::setw(4) << _pairs[i][0] <<
                    std::setw(4) << "-" <<
                    std::setw(4) << AtomType::lookupAtomSym(_atoms[_pairs[i][1] - 1]) << std::setw(4) << _pairs[i][1] << ":" <<
                    std::setprecision(5) << std::fixed <<
                    std::setw(13) << _interCoor->_Mdata[i + _nBonds + _nAngles + _nDihedrals] << 
                    std::setw(14) << _e->_Mdata[i + _nBonds + _nAngles + _nDihedrals] << std::endl;
        file << "Analytical gradient of overall energy: (in kcal/mol/Å)" << std::endl;
        for (size_t i = 0; i < _nAtoms; ++i)
            file << AtomType::lookupAtomSym(_atoms[i]) <<
                    std::setprecision(6) << std::fixed <<
                    std::setw(14) << _gx->_Mdata[i] <<
                    std::setw(14) << _gx->_Mdata[i + _nAtoms] <<
                    std::setw(14) << _gx->_Mdata[i + _nAtoms * 2] << std::endl;
        Matrix gxStretch(DIMENSION * _nAtoms, 1);
        double start = 0, end = _nBonds - 1;
        gxStretch = _B->cutR(start, end).transpose() * _gq->cutR(start, end);
        file << "Analytical gradient of stretching energy:" << std::endl;
        for (size_t i = 0; i < _nAtoms; ++i)
            file << AtomType::lookupAtomSym(_atoms[i]) <<
                    std::setprecision(6) << std::fixed <<
                    std::setw(14) << gxStretch._Mdata[i] <<
                    std::setw(14) << gxStretch._Mdata[i + _nAtoms] <<
                    std::setw(14) << gxStretch._Mdata[i + _nAtoms * 2] << std::endl;
        Matrix gxBend(DIMENSION * _nAtoms, 1);
        start = end + 1; end = end + _nAngles;
        gxBend = _B->cutR(start, end).transpose() * _gq->cutR(start, end);
        file << "Analytical gradient of bending energy:" << std::endl;
        for (size_t i = 0; i < _nAtoms; ++i)
            file << AtomType::lookupAtomSym(_atoms[i]) <<
                    std::setprecision(6) << std::fixed <<
                    std::setw(14) << gxBend._Mdata[i] <<
                    std::setw(14) << gxBend._Mdata[i + _nAtoms] <<
                    std::setw(14) << gxBend._Mdata[i + _nAtoms * 2] << std::endl;
        Matrix gxTorsion(DIMENSION * _nAtoms, 1);
        start = end + 1; end = end + _nDihedrals;
        gxTorsion = _B->cutR(start, end).transpose() * _gq->cutR(start, end);                    
        file << "Analytical gradient of torsional energy:" << std::endl;
        for (size_t i = 0; i < _nAtoms; ++i)
            file << AtomType::lookupAtomSym(_atoms[i]) <<
                    std::setprecision(6) << std::fixed <<
                    std::setw(14) << gxTorsion._Mdata[i] <<
                    std::setw(14) << gxTorsion._Mdata[i + _nAtoms] <<
                    std::setw(14) << gxTorsion._Mdata[i + _nAtoms * 2] << std::endl;
        Matrix gxVDW(DIMENSION * _nAtoms, 1);
        start = end + 1; end = end + _nPairs;
        gxVDW = _B->cutR(start, end).transpose() * _gq->cutR(start, end);  
        file << "Analytical gradient of VDW energy:" << std::endl;
        for (size_t i = 0; i < _nAtoms; ++i)
            file << AtomType::lookupAtomSym(_atoms[i]) <<
                    std::setprecision(6) << std::fixed <<
                    std::setw(14) << gxVDW._Mdata[i] <<
                    std::setw(14) << gxVDW._Mdata[i + _nAtoms] <<
                    std::setw(14) << gxVDW._Mdata[i + _nAtoms * 2] << std::endl;
        
        break;
    }
    case 2:
    {
        _M = new Matrix(DIMENSION * _nAtoms, DIMENSION * _nAtoms);
        Hessian::guessM(*_M);
        file << "Geometry optimization gradient RMS threshold:" << 
                std::setprecision(8) << std::fixed <<
                std::setw(15) << tol << std::endl;
        file << "Gradient of the energy for initial structure (in kcal/mol/Angstrom)" << std::endl;
        for (size_t i = 0; i < _nAtoms; ++i)
            file << AtomType::lookupAtomSym(_atoms[i]) << std::setprecision(6) << std::fixed <<
            std::setw(14) << _gx->_Mdata[i] <<
            std::setw(14) << _gx->_Mdata[i + _nAtoms] <<
            std::setw(14) << _gx->_Mdata[i + _nAtoms * 2] << std::endl;
        file << "Initial Guess for the Inverse Hessian" << std::endl;
        for (size_t i = 0; i < _M->_row; ++i)
        {
            file << "Row" << std::setw(5) << i + 1 << "  " << std::setprecision(6) << std::fixed;
            for (size_t j = 0; j < _M->_column; ++j)
            {
                file << std::setw(10) << (*_M)(i,j);
            }
            file << std::endl;
        }        
        evolve(file, 2, tol);
        
        file << std::endl << "##########################" << std::endl;
        file <<              "# Optimization converged #" << std::endl;
        file <<              "##########################" << std::endl << std::endl;
        file << "Final energy at mimimum:" << std::setprecision(8) << std::fixed << std::setw(16) << _e->sum() << " kcal/mol" << std::endl;
        file << "Final coordinates:" << std::endl;
        for (size_t i = 0; i < _nAtoms; ++i)
            file << std::setw(7) << AtomType::lookupAtomSym(_atoms[i]) << std::setprecision(6) << std::fixed <<
            std::setw(14) << _pos->_Mdata[i] <<
            std::setw(14) << _pos->_Mdata[i + _nAtoms] <<
            std::setw(14) << _pos->_Mdata[i + _nAtoms * 2] << std::endl;
        break;
    }
    case 3:
    {
        file << "Geometry optimization gradient RMS threshold:" << 
                std::setprecision(8) << std::fixed <<
                std::setw(15) << tol << std::endl;
        file << "Internal potential energy" << std::endl << 
                std::setw(12) << std::setprecision(6) << std::fixed << 
                _e->sum() << std::endl;
        evolve(file, 3, tol);
        file << std::endl << "##########################" << std::endl;
        file <<              "# Optimization converged #" << std::endl;
        file <<              "##########################" << std::endl << std::endl;
        file << "Final energy at mimimum:" << std::setprecision(8) << std::fixed << std::setw(16) << _e->sum() << " kcal/mol" << std::endl;
        file << "Final coordinates:" << std::endl;
        for (size_t i = 0; i < _nAtoms; ++i)
            file << std::setw(7) << AtomType::lookupAtomSym(_atoms[i]) << std::setprecision(6) << std::fixed <<
            std::setw(14) << _pos->_Mdata[i] <<
            std::setw(14) << _pos->_Mdata[i + _nAtoms] <<
            std::setw(14) << _pos->_Mdata[i + _nAtoms * 2] << std::endl;
        break;
    }
    case 4:
    {
        _M = new Matrix(DIMENSION * _nAtoms, DIMENSION * _nAtoms);
        Hessian::guessM(*_M);
        file << "Geometry optimization gradient RMS threshold:" << 
                std::setprecision(8) << std::fixed <<
                std::setw(15) << tol << std::endl;
        file << "Gradient of the energy for initial structure (in kcal/mol/Angstrom)" << std::endl;
        for (size_t i = 0; i < _nAtoms; ++i)
            file << AtomType::lookupAtomSym(_atoms[i]) << std::setprecision(6) << std::fixed <<
            std::setw(14) << _gx->_Mdata[i] <<
            std::setw(14) << _gx->_Mdata[i + _nAtoms] <<
            std::setw(14) << _gx->_Mdata[i + _nAtoms * 2] << std::endl;
        evolve(file, 4, tol);
        file << std::endl << "##########################" << std::endl;
        file <<              "# Optimization converged #" << std::endl;
        file <<              "##########################" << std::endl << std::endl;
        file << "Final energy at mimimum:" << std::setprecision(8) << std::fixed << std::setw(16) << _e->sum() << " kcal/mol" << std::endl;
        file << "Final coordinates:" << std::endl;
        for (size_t i = 0; i < _nAtoms; ++i)
            file << std::setw(7) << AtomType::lookupAtomSym(_atoms[i]) << std::setprecision(6) << std::fixed <<
            std::setw(14) << _pos->_Mdata[i] <<
            std::setw(14) << _pos->_Mdata[i + _nAtoms] <<
            std::setw(14) << _pos->_Mdata[i + _nAtoms * 2] << std::endl;
        break;
    }
    case 6:
    {       
        file << "Geometry optimization gradient RMS threshold:" << 
                std::setprecision(8) << std::fixed <<
                std::setw(15) << tol << std::endl;
        file << "Internal potential energy" << std::endl << 
                std::setw(12) << std::setprecision(6) << std::fixed << 
                _e->sum() << std::endl;
        evolve(file, 6, tol);
        file << std::endl << "##########################" << std::endl;
        file <<              "# Optimization converged #" << std::endl;
        file <<              "##########################" << std::endl << std::endl;
        file << "Final energy at mimimum:" << std::setprecision(8) << std::fixed << std::setw(16) << _e->sum() << " kcal/mol" << std::endl;
        file << "Final coordinates:" << std::endl;
        for (size_t i = 0; i < _nAtoms; ++i)
            file << std::setw(7) << AtomType::lookupAtomSym(_atoms[i]) << std::setprecision(6) << std::fixed <<
            std::setw(14) << _pos->_Mdata[i] <<
            std::setw(14) << _pos->_Mdata[i + _nAtoms] <<
            std::setw(14) << _pos->_Mdata[i + _nAtoms * 2] << std::endl;
        break;
    }
    default:
        std::cout << "No output format " << format << std::endl;
        exit(1);
}
}

Molecule* Molecule::readFromFile(const char* const fileName)
{
    FILE *file = fopen(fileName, "r");
    if (file == nullptr)
    {
        std::cout << "The input file does NOT exist!" << std::endl;
        exit(1);
    }

    const int buffersize = 256;
    char *buffer = new char[buffersize];

    if (fgets(buffer, buffersize * sizeof(char), file) == nullptr)
    {
        std::cout << "The input file is EMPTY!" << std::endl;
        exit(1);
    }

    size_t nAtoms, nBonds;

    //Get the number of atoms and number of bonds
    if (sscanf(buffer, " %zu %zu", &nAtoms, &nBonds) != 2)
    {
        std::cout << "The input file does NOT have the number of atoms and the number of bonds!" << std::endl;
        exit(1);
    }

    double *pos = new double[DIMENSION * nAtoms];

    unsigned int *atoms = new unsigned int[nAtoms];
    unsigned int **bonds = new unsigned int*[nBonds];
    for (size_t i = 0; i < nBonds; ++i) 
        bonds[i] = new unsigned int[N_ATOM_STRETCH];



    //Read the position and symbol of atoms
    for (size_t i = 0; i < nAtoms; ++i)
    {
        fgets(buffer, buffersize * sizeof(char), file);
        char atomType[2];
        sscanf(buffer, " %lf %lf %lf %s", pos + i, pos + nAtoms + i, pos + 2 * nAtoms + i, atomType);
        atoms[i] = AtomType::lookupAtomNumber(atomType);
    }

    for (size_t i = 0; i < nBonds; ++i)
    {
        fgets(buffer, buffersize * sizeof(char), file);
        sscanf(buffer, " %u %u", &bonds[i][0], &bonds[i][1]);
    }
    Matrix *p = new Matrix(DIMENSION * nAtoms, 1, pos);

    Molecule *mol = new Molecule(nAtoms, nBonds, atoms, bonds, p);


    delete[] buffer;
    fclose(file);
    return mol;
}


void Molecule::findForceConst()
{
    _para1 = new Matrix(_nBonds + _nAngles + _nDihedrals + _nPairs, 1);
    _para2 = new Matrix(_nBonds + _nAngles + _nDihedrals + _nPairs, 1);
    MolForcePara stretch(MolForcePara::fSTRETCH);
    for (size_t i = 0; i < _nBonds; ++i)
    {
        std::vector<unsigned int> ats;
        for (size_t j = 0; j < N_ATOM_STRETCH; ++j)
            ats.push_back(_atoms[_bonds[i][j] - 1]);
        std::array<double,2> p = stretch.findPara(ats);
        _para1->_Mdata[i] = p[0];
        _para2->_Mdata[i] = p[1];
    }
    MolForcePara bend(MolForcePara::fBEND);
    for (size_t i = 0; i < _nAngles; ++i)
    {
        std::vector<unsigned int> ats;
        for (size_t j = 0; j < N_ATOM_BEND; ++j)
            ats.push_back(_atoms[_angles[i][j] - 1]);
        std::array<double,2> p = bend.findPara(ats);
        _para1->_Mdata[i + _nBonds] = p[0];
        _para2->_Mdata[i + _nBonds] = p[1];
    }
    MolForcePara torsion(MolForcePara::fTORSION);
    for (size_t i = 0; i < _nDihedrals; ++i)
    {
        std::vector<unsigned int> ats;
        for (size_t j = 0; j < N_ATOM_TORSION; ++j)
            ats.push_back(_atoms[_dihedrals[i][j] - 1]);
        std::array<double,2> p = torsion.findPara(ats);
        _para1->_Mdata[i + _nBonds + _nAngles] = p[0];
        _para2->_Mdata[i + _nBonds + _nAngles] = p[1];
    }
    MolForcePara vdw(MolForcePara::fVDW);
    for (size_t i = 0; i < _nPairs; ++i)
    {
        // std::cout << _pairs[i][0] << _pairs[i][1] << std::endl;
        bool ifVDW = true;
        for (size_t j = 0; j < _nAngles; ++j)
        {
            if (findInArray(_angles[j], N_ATOM_BEND, _pairs[i][0])
             && findInArray(_angles[j], N_ATOM_BEND, _pairs[i][1]))
            {
                ifVDW = false;
                break;
            }
        }
        if (ifVDW)
        {
            std::vector<unsigned int> ats;
            for (size_t j = 0; j < N_ATOM_VDW; ++j)
                ats.push_back(_atoms[_pairs[i][j] - 1]);
            std::array<double,2> p = vdw.findPara(ats);
            _para1->_Mdata[i + _nBonds + _nAngles + _nDihedrals] = p[0];
            _para2->_Mdata[i + _nBonds + _nAngles + _nDihedrals] = p[1];
        }
        else
        {
            _para1->_Mdata[i + _nBonds + _nAngles + _nDihedrals] = 0.0;
            _para2->_Mdata[i + _nBonds + _nAngles + _nDihedrals] = 0.0;
        }
    }
}

void Molecule::calEnergy(Matrix &e, Matrix &interC)
{
    size_t start = 0, end = _nBonds - 1;
    MolForcePara stretch(MolForcePara::fSTRETCH);
    stretch.calEnergy(e.cutR(start, end), interC.cutR(start, end), _para1->cutR(start, end), _para2->cutR(start, end));
    MolForcePara bend(MolForcePara::fBEND);
    start = end + 1; end = end + _nAngles;
    bend.calEnergy(e.cutR(start, end), interC.cutR(start, end), _para1->cutR(start, end), _para2->cutR(start, end));
    MolForcePara torsion(MolForcePara::fTORSION);
    start = end + 1; end = end + _nDihedrals;
    torsion.calEnergy(e.cutR(start, end), interC.cutR(start, end), _para1->cutR(start, end), _para2->cutR(start, end));
    MolForcePara vdw(MolForcePara::fVDW);
    start = end + 1; end = end + _nPairs;
    vdw.calEnergy(e.cutR(start, end), interC.cutR(start, end), _para1->cutR(start, end), _para2->cutR(start, end));
    return;
}

void Molecule::calGradient(Matrix &gq, Matrix &interC)
{
    size_t start = 0, end = _nBonds - 1;
    MolForcePara stretch(MolForcePara::fSTRETCH);
    stretch.calGradient(gq.cutR(start, end), interC.cutR(start, end), _para1->cutR(start, end), _para2->cutR(start, end));
    MolForcePara bend(MolForcePara::fBEND);
    start = end + 1; end = end + _nAngles;
    bend.calGradient(gq.cutR(start, end), interC.cutR(start, end), _para1->cutR(start, end), _para2->cutR(start, end));
    MolForcePara torsion(MolForcePara::fTORSION);
    start = end + 1; end = end + _nDihedrals;
    torsion.calGradient(gq.cutR(start, end), interC.cutR(start, end), _para1->cutR(start, end), _para2->cutR(start, end));
    MolForcePara vdw(MolForcePara::fVDW);
    start = end + 1; end = end + _nPairs;
    vdw.calGradient(gq.cutR(start, end), interC.cutR(start, end), _para1->cutR(start, end), _para2->cutR(start, end));
    return;
}

void Molecule::calInterCoord(const Matrix &posA, Matrix &interC)
{
    Internal stretch(Internal::STRETCH);
    for(size_t i = 0; i < _nBonds; ++i)
    {
        std::vector<Matrix> pos;
        for (size_t j = 0; j < N_ATOM_STRETCH; ++j)
        {
            pos.push_back(Matrix(DIMENSION, 1));
            for (size_t k = 0; k < DIMENSION; ++k)
                pos[j]._Mdata[k] = posA._Mdata[_bonds[i][j] - 1 + k * _nAtoms];
        }
        interC._Mdata[i] = stretch.calculateIC(pos);
    }
    Internal bend(Internal::BEND);
    for(size_t i = 0; i < _nAngles; ++i)
    {
        std::vector<Matrix> pos;
        for (size_t j = 0; j < N_ATOM_BEND; ++j)
        {
            pos.push_back(Matrix(DIMENSION, 1));
            for (size_t k = 0; k < DIMENSION; ++k)
                pos[j]._Mdata[k] = posA._Mdata[_angles[i][j] - 1 + k * _nAtoms];
        }
        interC._Mdata[i + _nBonds] = bend.calculateIC(pos);
    }
    Internal torsion(Internal::TORSION);
    for(size_t i = 0; i < _nDihedrals; ++i)
    {
        std::vector<Matrix> pos;
        for (size_t j = 0; j < N_ATOM_TORSION; ++j)
        {
            pos.push_back(Matrix(DIMENSION, 1));
            for (size_t k = 0; k < DIMENSION; ++k)
                pos[j]._Mdata[k] = posA._Mdata[_dihedrals[i][j] - 1 + k * _nAtoms];
        }
        interC._Mdata[i + _nBonds + _nAngles] = torsion.calculateIC(pos);
    }
    Internal vdw(Internal::VDW);
    for(size_t i = 0; i < _nPairs; ++i)
    {
        std::vector<Matrix> pos;
        for (size_t j = 0; j < N_ATOM_VDW; ++j)
        {
            pos.push_back(Matrix(DIMENSION, 1));
            for (size_t k = 0; k < DIMENSION; ++k)
                pos[j]._Mdata[k] = posA._Mdata[_pairs[i][j] - 1 + k * _nAtoms];
        }
        interC._Mdata[i + _nBonds + _nAngles + _nDihedrals] = vdw.calculateIC(pos);
        std::vector<double> driv = vdw.calculateB(pos);
    }
}


void Molecule::calInterCoord(const Matrix &posA, Matrix &interC, Matrix &B)
{
    B.setZero();
    Internal stretch(Internal::STRETCH);
    for(size_t i = 0; i < _nBonds; ++i)
    {
        std::vector<Matrix> pos;
        for (size_t j = 0; j < N_ATOM_STRETCH; ++j)
        {
            pos.push_back(Matrix(DIMENSION, 1));
            for (size_t k = 0; k < DIMENSION; ++k)
                pos[j]._Mdata[k] = posA._Mdata[_bonds[i][j] - 1 + k * _nAtoms];
        }
        interC._Mdata[i] = stretch.calculateIC(pos);
        std::vector<double> driv = stretch.calculateB(pos);
        size_t count = 0;
        for (size_t j = 0; j < N_ATOM_STRETCH; ++j)
        {
            for (size_t k = 0; k < DIMENSION; ++k)
            {
                B(i, _bonds[i][j] - 1 + k * _nAtoms) = driv[count];
                ++count;
            }
        }
    }
    Internal bend(Internal::BEND);
    for(size_t i = 0; i < _nAngles; ++i)
    {
        std::vector<Matrix> pos;
        for (size_t j = 0; j < N_ATOM_BEND; ++j)
        {
            pos.push_back(Matrix(DIMENSION, 1));
            for (size_t k = 0; k < DIMENSION; ++k)
                pos[j]._Mdata[k] = posA._Mdata[_angles[i][j] - 1 + k * _nAtoms];
        }
        interC._Mdata[i + _nBonds] = bend.calculateIC(pos);
        std::vector<double> driv = bend.calculateB(pos);
        size_t count = 0;
        for (size_t j = 0; j < N_ATOM_BEND; ++j)
        {
            for (size_t k = 0; k < DIMENSION; ++k)
            {
                B(i + _nBonds, _angles[i][j] - 1 + k * _nAtoms) = driv[count];
                ++count;
            }
        }
    }
    Internal torsion(Internal::TORSION);
    for(size_t i = 0; i < _nDihedrals; ++i)
    {
        std::vector<Matrix> pos;
        for (size_t j = 0; j < N_ATOM_TORSION; ++j)
        {
            pos.push_back(Matrix(DIMENSION, 1));
            for (size_t k = 0; k < DIMENSION; ++k)
                pos[j]._Mdata[k] = posA._Mdata[_dihedrals[i][j] - 1 + k * _nAtoms];
        }
        interC._Mdata[i + _nBonds + _nAngles] = torsion.calculateIC(pos);
        std::vector<double> driv = torsion.calculateB(pos);
        size_t count = 0;
        for (size_t j = 0; j < N_ATOM_TORSION; ++j)
        {
            for (size_t k = 0; k < DIMENSION; ++k)
            {
                B(i + _nBonds + _nAngles, _dihedrals[i][j] - 1 + k * _nAtoms) = driv[count];
                ++count;
            }
        }
    }
    Internal vdw(Internal::VDW);
    for(size_t i = 0; i < _nPairs; ++i)
    {
        std::vector<Matrix> pos;
        for (size_t j = 0; j < N_ATOM_VDW; ++j)
        {
            pos.push_back(Matrix(DIMENSION, 1));
            for (size_t k = 0; k < DIMENSION; ++k)
                pos[j]._Mdata[k] = posA._Mdata[_pairs[i][j] - 1 + k * _nAtoms];
        }
        interC._Mdata[i + _nBonds + _nAngles + _nDihedrals] = vdw.calculateIC(pos);
        std::vector<double> driv = vdw.calculateB(pos);
        size_t count = 0;
        for (size_t j = 0; j < N_ATOM_VDW; ++j)
        {
            for (size_t k = 0; k < DIMENSION; ++k)
            {
                B(i + _nBonds + _nAngles + _nDihedrals, _pairs[i][j] - 1 + k * _nAtoms) = driv[count];
                ++count;
            }
        }
    }
}


void Molecule::evolve(std::ofstream &file, unsigned int mode, double tol)
{
    switch(mode) {
    /*
        1 : singlepoint energy ( so it do not enter this function )
        2 : cartesian optimization with debugging output
        3 : internal optimization with debugging output
        4 : cartesian optimization with less output
        6 : internal optimization with less output
    */
    case 2:
    {
        double RMS, a;
        Matrix y(_gx->_row, _gx->_column);
        Matrix s(_gx->_row, _gx->_column);
        Matrix v(_gx->_row, _gx->_column);
        unsigned int optCount = 0;
        do
        {
            ++optCount;
            Matrix p = - (*_M) * (*_gx);
            file << std::endl << "***** Geometry optimization cycle number" << std::setw(4) << optCount << " ******" << std::endl << std::endl;
            file << "Predicted structure change p_k = -M_k grad V(r_k)" << std::endl;
            for (size_t i = 0; i < _nAtoms; ++i)
                file << AtomType::lookupAtomSym(_atoms[i]) << std::setprecision(6) << std::fixed <<
                        std::setw(14) << p._Mdata[i] <<
                        std::setw(14) << p._Mdata[i + _nAtoms] <<
                        std::setw(14) << p._Mdata[i + _nAtoms * 2] << std::endl;
            a = lineSearch(file, p, 1.0, true);
            s = a * p;
            calGradient(*_gq, *_interCoor);
            Matrix gxOld = *_gx;
            *_gx = _B->transpose() * (*_gq);
            y = *_gx - gxOld;
            Hessian::updateM(*_M, s, y);
            RMS = norm(*_gx) / pow(_nAtoms * DIMENSION, 0.5);
            file << " And GRMS:" << std::setprecision(7) << std::fixed << std::setw(15) << RMS << std::endl;
            file << "New gradient g_k+1:" << std::endl;
            for (size_t i = 0; i < _nAtoms; ++i)
                file << AtomType::lookupAtomSym(_atoms[i]) << std::setprecision(6) << std::fixed <<
                std::setw(14) << _gx->_Mdata[i] <<
                std::setw(14) << _gx->_Mdata[i + _nAtoms] <<
                std::setw(14) << _gx->_Mdata[i + _nAtoms * 2] << std::endl;
            file << "Updated guess for the Inverse Hessian (equation 12)" << std::endl;
            for (size_t i = 0; i < _M->_row; ++i) {
                file << "Row" << std::setw(5) << i + 1 << "  " << std::setprecision(6) << std::fixed;
                for (size_t j = 0; j < _M->_column; ++j) {
                    file << std::setw(10) << (*_M)(i,j);
                }
                file << std::endl;
            }
        } while(RMS > tol);
        break;
    }
    case 3:
    {
        double RMS;
        unsigned int optCount = 0;
        Matrix B = _B->cutR(0, _nBonds + _nAngles + _nDihedrals - 1);
        file << "Wilson B Matrix at the initial structure:" << std::setw(4) << B._row << " rows of" << std::setw(4) << B._column << " elements" << std::endl;
        for (size_t i = 0; i < B._row; ++i)
        {
            file << std::setw(4) << i + 1 << ":";
            for (size_t j = 0; j < B._column; ++j) file << std::setprecision(5) << std::fixed << std::setw(12) << B(i, j);
            file << std::endl;
        }
        Matrix G = B * B.transpose();
        file << "G matrix, the product of B with its own transpose (square, dimension" << std::setw(4) << G._row << ")" << std::endl;
        for (size_t i = 0; i < G._row; ++i)
        {
            file << std::setw(4) << i + 1 << ":";
            for (size_t j = 0; j < G._column; ++j) file << std::setprecision(5) << std::fixed << std::setw(12) << G(i, j);
            file << std::endl;
        }
        Matrix invG = G.getInv();

        file << "Inverse G Matrix (square matrix with dimension" << std::setw(4) << invG._row << ") at the initial structure:" << std::endl;
        for (size_t i = 0; i < invG._row; ++i)
        {
            file << std::setw(4) << i + 1 << ":";
            for (size_t j = 0; j < invG._column; ++j) file << std::setprecision(5) << std::fixed << std::setw(12) << invG(i, j);
            file << std::endl;
        }

        file << "Initial gradient in terms of the internal coordinates (kcal/mol/Angstrom or kcal/mol/radian):" << std::endl;

        Matrix gq = invG * B * (*_gx);
        for (size_t i = 0; i < gq._row; ++i) file << std::setprecision(6) << std::fixed << std::setw(12) << gq(i,0);
        file << std::endl;

        _M = new Matrix(gq._row, gq._row);
        Hessian::guessM(*_M, _nBonds, _nAngles, _nDihedrals);
        file << "Initial guess for the inverse Hessian M in internal coordinates:" << std::endl;
        for (size_t i = 0; i < invG._row; ++i)
        {
            file << std::setw(4) << i + 1 << ":";
            for (size_t j = 0; j < invG._column; ++j) file << std::setprecision(5) << std::fixed << std::setw(12) << (*_M)(i, j);
            file << std::endl;
        }
        
        file << std::endl << "##################################" << std::endl;
        file <<              "# Start of Geometry Optimization #" << std::endl;
        file <<              "##################################" << std::endl << std::endl;

        const double stepSize = 0.020;
        file << "Reminder: maximum step size is" << std::setw(8) << stepSize << std::endl << std::endl;
        do
        {
            ++optCount;
            file << "*** Optimization Cycle" << std::setw(4) << optCount << " ***" << std::endl << std::endl;

            Matrix p = - (*_M) * gq;
            file << "Predicted update step in internal coordinates s_k (prior to possible scaling):" << std::endl;
            for (size_t i = 0; i < p._row; ++i) file << std::setprecision(6) << std::fixed << std::setw(12) << p(i,0);
            file << std::endl;
            double pRMS = norm(p) / pow(p._row, 0.5);
            if (pRMS > stepSize)
            {
                file << "Predicted step is too long: RMS length:" << std::setw(13) << pRMS << std::endl;
                p = (stepSize / pRMS) * p;
                file << "Scaled update step in internal coordinates s_k:" << std::endl;
                for (size_t i = 0; i < p._row; ++i) file << std::setprecision(6) << std::fixed << std::setw(12) << p(i,0);
                file << std::endl;
            }
            Matrix q = _interCoor->cutR(0, _nBonds + _nAngles + _nDihedrals - 1);
            Matrix qOld = q;
            q = q + p;
            Matrix x = *_pos;
            Matrix q1 = *_interCoor;
            x = x + B.transpose() * invG * p;
            file << "Initially predicted dx = BTG-sk:" << std::endl;
            for (size_t i = 0; i < x._row; ++i) file << std::setprecision(6) << std::fixed << std::setw(12) << x(i,0) - (*_pos)(i,0);
            file << std::endl;
            unsigned int cartCount = 0;
            double dxThresh = 1e-5;
            double maxDx;
            do
            {
                Matrix oldX = x;
                ++cartCount;
                calInterCoord(x,q1);
                file << "  Cartesian fitting iteration number" << std::setw(4) << cartCount << std::endl;
                file << "current set of internals q_(k+1)^(j+1):" << std::endl;
                for (size_t i = 0; i < _nBonds + _nAngles + _nDihedrals; ++i) file << std::setprecision(6) << std::fixed << std::setw(12) << q1(i,0);
                file << std::endl;
                p = q - q1.cutR(0, _nBonds + _nAngles + _nDihedrals - 1);
                for (size_t i = _nBonds + _nAngles; i < _nBonds + _nAngles + _nDihedrals; ++i) p(i,0) = fabs(p(i,0)) > PI ? p(i,0) - ((p(i,0) > 0) - (p(i,0) < 0)) * 2 * PI : p(i, 0);
                file << "difference between these internals and the desired internals, s_q,k^j+1:" << std::endl;
                for (size_t i = 0; i < p._row; ++i) file << std::setprecision(6) << std::fixed << std::setw(12) << p(i,0);
                file << std::endl;
                Matrix dx = B.transpose() * invG * p;
                x = x + dx;
                file << "Corresponding Cartesians x+k+1^j:" << std::endl;
                for (size_t i = 0; i < x._row; ++i) file << std::setprecision(6) << std::fixed << std::setw(12) << x(i,0);
                file << std::endl;
                maxDx = 0;
                for (size_t i = 0; i < dx._row; ++i) if (fabs(dx(i,0)) > maxDx) maxDx = fabs(dx(i,0)); 
                file << "Maximum change in x from previous iteration" << std::setw(13) << maxDx << std::endl;
            } while (maxDx > dxThresh);
            file << "  fitting of Cartesian coordinates considered converged." << std::endl << std::endl;
            *_pos = x;
            calInterCoord(*_pos, *_interCoor, *_B);
            calGradient(*_gq, *_interCoor);

            file << "New coordiantes:" << std::endl;
            for (size_t i = 0; i < _nAtoms; ++i)
                file << "      " << AtomType::lookupAtomSym(_atoms[i]) << std::setprecision(6) << std::fixed <<
                std::setw(14) << _pos->_Mdata[i] <<
                std::setw(14) << _pos->_Mdata[i + _nAtoms] <<
                std::setw(14) << _pos->_Mdata[i + _nAtoms * 2] << std::endl;
            B = _B->cutR(0, _nBonds + _nAngles + _nDihedrals - 1);
            G = B * B.transpose();
            invG = G.getInv();              
            file << "New set of internals q (note these are the ones that correspond to the best fit Cartesians):" << std::endl;
                for (size_t i = 0; i < _nBonds + _nAngles + _nDihedrals; ++i) file << std::setprecision(6) << std::fixed << std::setw(12) << (*_interCoor)(i,0);
            file << std::endl;
            file << "Wilson B Matrix at the new structure:" << std::endl;
            for (size_t i = 0; i < B._row; ++i)
            {
                file << std::setw(4) << i + 1 << ":";
                for (size_t j = 0; j < B._column; ++j) file << std::setprecision(5) << std::fixed << std::setw(12) << B(i, j);
                file << std::endl;
            }
            file << "Inverse G Matrix at the new structure:" << std::endl;
            for (size_t i = 0; i < invG._row; ++i)
            {
                file << std::setw(4) << i + 1 << ":";
                for (size_t j = 0; j < invG._column; ++j) file << std::setprecision(5) << std::fixed << std::setw(12) << invG(i, j);
                file << std::endl;
            }
            (*_gx) = _B->transpose() * (*_gq);
            Matrix gqOld = gq;
            q = _interCoor->cutR(0, _nBonds + _nAngles + _nDihedrals - 1);
            gq = invG * B * (*_gx);
            RMS = norm(*_gx) / std::pow(_gx->_row, 0.5);
            file << "Gradient in terms of the internal coordinates:" << std::endl;
            for (size_t i = 0; i < gq._row; ++i) file << std::setprecision(6) << std::fixed << std::setw(12) << gq(i,0);
            file << std::endl;
            file << "Old and new energies:" << std::setw(13) << _e->sum();
            calEnergy(*_e, *_interCoor);
            file << std::setw(13) << _e->sum() << " And GRMS:" << std::setw(13) << RMS << std::endl << std::endl;
            if (RMS < tol) break;
            Matrix y = gq - gqOld;
            Matrix v = (*_M) * y;
            Matrix s = q - qOld;
            for (size_t i = _nBonds + _nAngles; i < _nBonds + _nAngles + _nDihedrals; ++i) s(i,0) = fabs(s(i,0)) > PI ? s(i,0) - ((s(i,0) > 0) - (s(i,0) < 0)) * 2 * PI : s(i, 0);
            file << "Change in the gradient in terms of the internal coordinates: (vector yk)" << std::endl;
            y.transpose().print(file);
            file << "Change in the structure in terms of the internal coordinates: (vector sk)" << std::endl;
            s.transpose().print(file);
            file << "Product of the inverse Hessian and the change in gradient, (vector vk):" << std::endl;
            v.transpose().print(file);
            file << "Dot product of sk with yk:" << std::endl;
            double dsy = dot(s, y), dyv = dot(y, v);
            file << std::setprecision(12) << std::fixed << dot(s, y) << std::endl;
            file << "Dot product of yk with vk:" << std::endl;
            file << std::setprecision(12) << std::fixed << dot(y, v) << std::endl;
            Matrix r1 = ((dsy + dyv) / std::pow(dsy, 2)) * s * s.transpose();
            file << "First term on right-hand side of eq. (20), (s∙y + y∙v)/((s∙y)^2) × s ⊗ s :" << std::endl;
            r1.print(file);
            Matrix r2 = (v * s.transpose() + s * v.transpose()) / dsy ;
            file << "Second term on right-hand side of eq. (20), (v ⊗ s + s ⊗ v) / (s∙y) :" << std::endl;
            r2.print(file);
            Hessian::updateM(*_M, s, y);
            file << "Updated approximate inverse Hessian M in internal coordinates:" << std::endl;
            for (size_t i = 0; i < _M->_row; ++i)
            {
                file << std::setw(4) << i + 1 << ":";
                for (size_t j = 0; j < _M->_column; ++j) file << std::setprecision(6) << std::fixed << std::setw(12) << (*_M)(i, j);
                file << std::endl;
            }

        } while (RMS > tol);
        break;
    }
    case 4:
    {
        double RMS, a;
        Matrix y(_gx->_row, _gx->_column);
        Matrix s(_gx->_row, _gx->_column);
        Matrix v(_gx->_row, _gx->_column);
        unsigned int optCount = 0;
        do
        {
            ++optCount;
            Matrix p = - (*_M) * (*_gx);
            file << std::endl << "***** Geometry optimization cycle number" << std::setw(4) << optCount << " ******" << std::endl << std::endl;
            a = lineSearch(file, p);
            s = a * p;
            calGradient(*_gq, *_interCoor);
            Matrix gxOld = *_gx;
            *_gx = _B->transpose() * (*_gq);
            y = *_gx - gxOld;
            Hessian::updateM(*_M, s, y);
            RMS = norm(*_gx) / pow(_nAtoms * DIMENSION, 0.5);
            file << " And GRMS:" << std::setprecision(7) << std::fixed << std::setw(15) << RMS << std::endl;
        } while(RMS > tol);
        break;
    }
    case 6:
    {
        double RMS;
        unsigned int optCount = 0;
        Matrix B = _B->cutR(0, _nBonds + _nAngles + _nDihedrals - 1);
        Matrix G = B * B.transpose();
        // inverse wilson B matrix
        Matrix invG = G.getInv();
        file << "Initial gradient in terms of the internal coordinates (kcal/mol/Angstrom or kcal/mol/radian):" << std::endl;

        // calculate gradient w.r.t internal coordinates
        Matrix gq = invG * B * (*_gx);
        for (size_t i = 0; i < gq._row; ++i) file << std::setprecision(6) << std::fixed << std::setw(12) << gq(i,0);
        file << std::endl;

        _M = new Matrix(gq._row, gq._row);
        Hessian::guessM(*_M, _nBonds, _nAngles, _nDihedrals);
        file << std::endl << "##################################" << std::endl;
        file <<              "# Start of Geometry Optimization #" << std::endl;
        file <<              "##################################" << std::endl << std::endl;
        const double stepSize = 0.020;
        do
        {
            ++optCount;
            file << "*** Optimization Cycle" << std::setw(4) << optCount << " ***" << std::endl << std::endl;
            // predict coordinate change in internal coordinates
            Matrix p = - (*_M) * gq;
            file << "Predicted update step in internal coordinates s_k (prior to possible scaling):" << std::endl;
            for (size_t i = 0; i < p._row; ++i) file << std::setprecision(6) << std::fixed << std::setw(12) << p(i,0);
            file << std::endl;
            double pRMS = norm(p) / pow(p._row, 0.5);
            if (pRMS > stepSize) {
                p = (stepSize / pRMS) * p;
            }
            Matrix q = _interCoor->cutR(0, _nBonds + _nAngles + _nDihedrals - 1);
            Matrix qOld = q;
            q = q + p;
            Matrix x = *_pos;
            Matrix q1 = *_interCoor;
            x = x + B.transpose() * invG * p;
            unsigned int cartCount = 0;
            double dxThresh = 1e-5;
            double maxDx;
            // find the optimized cartesian change
            do
            {
                Matrix oldX = x;
                ++cartCount;
                calInterCoord(x,q1);
                p = q - q1.cutR(0, _nBonds + _nAngles + _nDihedrals - 1);
                for (size_t i = _nBonds + _nAngles; i < _nBonds + _nAngles + _nDihedrals; ++i) p(i,0) = fabs(p(i,0)) > PI ? p(i,0) - ((p(i,0) > 0) - (p(i,0) < 0)) * 2 * PI : p(i, 0);
                Matrix dx = B.transpose() * invG * p;
                x = x + dx;
                file << std::endl;
                maxDx = 0;
                for (size_t i = 0; i < dx._row; ++i) if (fabs(dx(i,0)) > maxDx) maxDx = fabs(dx(i,0)); 
            } while (maxDx > dxThresh);
            *_pos = x;

            // calculate new internal coordinates and gradient and energy
            calInterCoord(*_pos, *_interCoor, *_B);
            calEnergy(*_e, *_interCoor);
            calGradient(*_gq, *_interCoor);
            file << "New coordiantes:" << std::endl;
            for (size_t i = 0; i < _nAtoms; ++i)
                file << "      " << AtomType::lookupAtomSym(_atoms[i]) << std::setprecision(6) << std::fixed <<
                std::setw(14) << _pos->_Mdata[i] <<
                std::setw(14) << _pos->_Mdata[i + _nAtoms] <<
                std::setw(14) << _pos->_Mdata[i + _nAtoms * 2] << std::endl;
            B = _B->cutR(0, _nBonds + _nAngles + _nDihedrals - 1);
            G = B * B.transpose();
            invG = G.getInv();
            (*_gx) = _B->transpose() * (*_gq);
            Matrix gqOld = gq;
            q = _interCoor->cutR(0, _nBonds + _nAngles + _nDihedrals - 1);
            gq = invG * B * (*_gx);
            RMS = norm(*_gx) / std::pow(_gx->_row, 0.5);
            file << std::endl;
            file << "Old and new energies:" << std::setw(13) << _e->sum();
            file << std::setw(13) << _e->sum() << " And GRMS:" << std::setw(13) << RMS << std::endl << std::endl;
            if (RMS < tol) break;
            Matrix y = gq - gqOld;
            Matrix s = q - qOld;
            for (size_t i = _nBonds + _nAngles; i < _nBonds + _nAngles + _nDihedrals; ++i) s(i,0) = fabs(s(i,0)) > PI ? s(i,0) - ((s(i,0) > 0) - (s(i,0) < 0)) * 2 * PI : s(i, 0);
            Hessian::updateM(*_M, s, y);
        } while (RMS > tol);
    break;
    }
    default:
        std::cout << "No mode " << mode << std::endl;
        exit(1);
    }
}

double Molecule::lineSearch(std::ofstream &file, const Matrix& p, double a, bool ifOutput)
{
    const double c1 = 0.1;
    Matrix posT(_pos->_row, _pos->_column);
    Matrix eT(_e->_row, _e->_column);
    Matrix interCT(_interCoor->_row, _interCoor->_column);
    Matrix BT(_B->_row, _B->_column);
    double energyOld = _e->sum();
    double energyNew;
    double rightSide = c1 * dot(p, *_gx);
    do
    {
        a *= 0.8;
        posT = (*_pos) + a * p;
        calInterCoord(posT, interCT, BT);
        calEnergy(eT, interCT);
        energyNew = eT.sum();
        if (ifOutput)
            file << "    Line search: alpha and energy:" << 
                    std::setprecision(6) << std::fixed <<
                    std::setw(13) << a << 
                    std::setw(13) << energyNew << std::endl;
    } while(energyNew > (energyOld + a * rightSide));
    *_pos = posT;
    *_e = eT;
    *_B = BT;
    *_interCoor = interCT;
    file << "New structure r_k+1 = r_k + s_k = r_k + alpha_k p_k" << std::endl;
    for (size_t i = 0; i < _nAtoms; ++i)
        file << AtomType::lookupAtomSym(_atoms[i]) << std::setprecision(6) << std::fixed <<
        std::setw(14) << _pos->_Mdata[i] <<
        std::setw(14) << _pos->_Mdata[i + _nAtoms] <<
        std::setw(14) << _pos->_Mdata[i + _nAtoms * 2] << std::endl;
    file << "Energy before update and after:" << std::setprecision(8) << std::fixed << 
            std::setw(15) << energyOld << std::setw(15) << energyNew;
    return a;
}
