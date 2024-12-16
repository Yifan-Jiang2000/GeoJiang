#include "mathj.h"
#include <iostream>
#include <fstream>
#include "molecule.h"
#include <cmath>
#include <fstream>
#include <cstring>
#include <unordered_map>
#include <cstdio>
#include <cassert>
#include <cstddef>
#include <ctime>
#include "f2c.h"
#include "clapack.h"


std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(' ');
    size_t last = str.find_last_not_of(' ');
    if (str.back() == '\n' || str.back() == 'r') --last;
    return str.substr(first, (last - first + 1));
}


int main()
{

    std::unordered_map<std::string, std::string> keywords;
    FILE *set = fopen("setting.ini", "r");
    assert(set && "setting.ini does not exist");
    char line[256];
    // size_t bytes;
    while (fgets(line, sizeof(line), set))
    {
        if (line[0] == '\n' || line[0] == '#') continue;
        char *splitter = strchr(line,'=');
        if (splitter)
        {
            *splitter = '\0';
            std::string key = trim(line);
            std::string value = trim(splitter + 1);
            keywords[key] = value;
        }
    }

    int task, outputLevel;
    double tol;


    // std::cout << keywords["task"] << std::endl;
    if (keywords["task"] == "cartesian") task = 2;
    else if (keywords["task"] == "internal") task = 3;
    else if (keywords["task"] == "singlepoint") task = 1;
    if (keywords["output"] == "debug") outputLevel = 1;
    else if (keywords["output"] == "info") outputLevel = 2;

    if (keywords["eigenSolver"] == "jiang") eigenSolver = &Matrix::eigenDecomposition;
    else eigenSolver = &Matrix::eigenDecompositionLAPACK;

    if (keywords.find("tolerance") == keywords.end()) tol = 0.001;
    else tol = std::stod(keywords["tolerance"]);


    if (task != 1) task *= outputLevel;


    fclose(set);

    Molecule *mol = Molecule::readFromFile(keywords["inputFile"].c_str());
    std::ofstream file(keywords["outputFile"].c_str());
    mol->output(file, task, tol);
    file.close();

    return 1;
}


