#pragma once
#include <iostream>
#include <cstdio>
using namespace std;

class Input_d{
    public:
        Input_d(const char* FilePath);
        ~Input_d();
        double cutoff;
        double dr;
        int mesh;
        double* f;
};