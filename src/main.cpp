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


int main(int argc, char *argv[])
{
    // double temp[16] = {1, 2, 3, 4, 2, 0, 5, 6, 3, 5, 0, 7, 4, 6, 7, 0};
    // double temp[16] = {4, 1, 2, 3, 1, 3, 1, 2, 2, 1, 2, 1, 3, 2, 1, 1};
    // double temp[16] = {2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2};
    // double temp[16] = {2, 1, 1, 6, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2};
    Matrix sym(10,10);
    for (size_t i = 0; i < sym._row; ++i) {
        for (size_t j = i; j < sym._column; ++j) {
            sym(i, j) = rand() / 2147483647.0;
            sym(j, i) = sym(i, j);
        }
    }
    sym.print();
    long i = 10;
    long nb = 5;
    double E[9];
    double TAU[9];
    double W[200];
    long LDW = 40;
    dlatrd_("U", &i, &nb, sym._Mdata, &i, E, TAU, W, &LDW);
    sym.print();

    // Matrix a = Matrix(500,500);
    // Matrix b = Matrix(500,500);
    // for (size_t i = 0; i < a._row * a._column; ++i)
    // {
    //     a._Mdata[i] = rand();
    //     b._Mdata[i] = rand();
    // }
    // Matrix c = Matrix(500,500);
    
    // std::clock_t start = clock();
    // c = a * b;
    // // c.print();
    // // long i = 500;
    // // double alpha = 1.;
    // // double beta = 0.;
    // // dgemm_("N", "N", &i, &i, &i, &alpha, a._Mdata, &i, b._Mdata, &i, &beta, c._Mdata, &i);
    
    // std::clock_t end = clock();
    // std::cout << (double)(end - start) / CLOCKS_PER_SEC << std::endl;


    // Matrix m = Matrix(4, 4);

    // for (int i = 0; i < 4; ++i)
    // {
    //     for (int j = 0; j < 4; ++j)
    //     {
    //         m._Mdata[i * m._column + j] = temp[i * m._column + j];
    //     }
    // }

    // Matrix hessen(m._row,m._column);
    // Matrix trans(m._row,m._column);
    // m.getUpHessenberg(trans,hessen);
    // m.print();
    // hessen.print();
    // std::cout << m._Mdata << "  " << hessen._Mdata << std::endl;
    // Matrix Q(m._row,m._column);
    // Matrix R(m._row,m._column);

    // Matrix lambda(m._row,m._column);
    // Matrix v(m._row,m._column);
    // Matrix invM = m.getInv();
    // m.eigenDecomposition(lambda,v);
    // m.print();
    // lambda.print();
    // v.print();

    // (v.transpose() * m * v).print();
    // Matrix invLambda = lambda.getInvDiag();
    // Matrix invM = v * invLambda * v.transpose();
    // (invM * m).print();
    

    // std::cout << m.det() << std::endl;

    // m.QRDecomposition(Q, R);
    // m.print();
    // Q.print();
    // R.print();

    // Matrix a(m);
    // a = m.guassEliminaton();
    // m.print();
    // a.print();
    // Matrix eigenVec(m._row, 1);
    // double lambda = m.powerIter(eigenVec);
    // m.print();
    // std::cout << m.det() << std::endl;
    // eigenVec.print();
    // std::cout << lambda << std::endl;

    // m.print();
    // std::cout << m.det() << std::endl;
    // m.transpose().print();
    // m.inverse().print();
    // return 1;





    // std::unordered_map<std::string, std::string> keywords;
    // FILE *set = fopen("setting.ini", "r");
    // assert(set && "setting.ini does not exist");
    // char line[256];
    // size_t bytes;
    // while (fgets(line, sizeof(line), set))
    // {
    //     if (line[0] == '\n' || line[0] == '#') continue;
    //     char *splitter = strchr(line,'=');
    //     if (splitter)
    //     {
    //         *splitter = '\0';
    //         std::string key = trim(line);
    //         std::string value = trim(splitter + 1);
    //         keywords[key] = value;
    //     }
    // }

    // int task, outputLevel;


    // // std::cout << keywords["task"] << std::endl;
    // if (keywords["task"] == "cartesian") task = 2;
    // else if (keywords["task"] == "internal") task = 3;
    // else if (keywords["task"] == "singlepoint") task = 1;
    // if (keywords["output"] == "debug") outputLevel = 1;
    // else if (keywords["output"] == "info") outputLevel = 2;



    // if (task != 1) task *= outputLevel;


    // fclose(set);

    // Molecule *mol = Molecule::readFromFile(keywords["inputFile"].c_str());
    // std::ofstream file(keywords["outputFile"].c_str());
    // mol->output(file, task);
    // file.close();

    // return 1;
}


