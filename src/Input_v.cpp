#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
//#include <omp.h>
#include <algorithm>
#include "Input_v.h"
using namespace std;

Input_v::Input_v(const char* FilePath){
    ifstream InputFile(FilePath);;
    /*stringstream InputStream;
    try{
        InputFile.open(FilePath);
        InputStream << InputFile.rdbuf();
        InputFile.close();
    }
    catch (ifstream::failure& e) { 
        cout << "ERROR::FILE_NOT_SUCCESSFULLY_READ:" << e.what() << endl;
    }*/
    string ASentence;
    while (getline(InputFile,ASentence)){
        int pos1 = ASentence.find(" ");
        string Premeter = ASentence.substr(0, pos1);
        int pos2 = ASentence.find(" ", pos1 + 1);
        string Operation = ASentence.substr(pos1 + 1, pos2 - pos1);
        if (Premeter == "nx") nx = stoi(Operation);
        if (Premeter == "ny") ny = stoi(Operation);
        if (Premeter == "nz") nz = stoi(Operation);
        if (Premeter == "V:"){
            //int nthreads = 4;
            //omp_set_num_threads(nthreads);
            venergy = new double** [nx];
            //#pragma omp parallel
            //{
            //#pragma omp for
            for (int i = 0 ; i < nx ; i ++){
                //#pragma omp critical
                venergy[i] = new double* [ny];
            }
            //#pragma omp for    
            for (int i = 0 ; i < nx ; i ++){
                for (int j = 0 ; j < ny ; j ++)
                    //#pragma omp critical
                    venergy[i][j] = new double [nz];
            }
            for (int i = 0 ; i < nx ; i ++)
                for (int j = 0 ; j < ny ; j ++)
                    for (int k = 0 ; k < nz ; k ++)
                        InputFile >> venergy[i][j][k];
                         //fscanf(InputFile, "%lf", &venergy[i][j][k]);
            //}
            
        }
        ngrid = nx * ny * nz;
    }
}

Input_v::~Input_v(){
    for (int i = 0 ; i < nx ; i ++)
        for (int j = 0 ; j < ny ; j ++){
            delete[] venergy[i][j];
            venergy[i][j] = NULL;
        }
    for (int i = 0 ; i < nx ; i ++){
        delete[] venergy[i];
        venergy[i] = NULL;
    }
    delete[] venergy;
    venergy = NULL;
}