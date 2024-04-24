#pragma once
#include <iostream>
#include <cstdio>
using namespace std;

class Input_p{
    public:
        Input_p(const char* FilePath);
        ~Input_p();
        int n;
        double** pos;
};