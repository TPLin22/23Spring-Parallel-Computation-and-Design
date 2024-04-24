#pragma once
#include <iostream>
#include <cstdio>
using namespace std;

class mysca{
    public:
        mysca(double ** H,int M,int myrank_mpi,int nprocs_mpi);
        ~mysca();
        double *A;
        double *Z;
        int* mA_array;
        int* nA_array;
        int* myrow_array;
        int* mycol_array;
        int *hasm,*hasn;
        double *w;
        double *work;
        double *GZ;
};