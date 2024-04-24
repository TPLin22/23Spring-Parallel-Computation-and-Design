#pragma once
#include <iostream>
#include <cstdio>
using namespace std;

class Input_v{
    public:
        Input_v(const char* FilePath);
        ~Input_v();
        int nx,ny,nz,ngrid;
        double*** venergy;
};