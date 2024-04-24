#pragma once
#include <cstdio>
#include <string>

using namespace std;

class Input_0{
    public:
        Input_0(const char* FilePath);
        ~Input_0();
        string diago_lib,points_path,v_path,distribution_path;
        double lx,ly,lz;
};
