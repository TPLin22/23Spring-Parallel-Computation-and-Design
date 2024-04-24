#pragma once
#include <iostream>
#include <cstdio>
using namespace std;

class mathwork{
    public:
        mathwork(double *f, int size_of_f, double dr);
        ~mathwork();
        double fl(double x);
        double calculat(double x,double y,double z,double x1[],double x2[],double v);
        double *xdata,*ydata;
        int size;
        double ** a,*h,*aver,*lan,*d,*nu,*beita,*yy,*xx;
        void uselapack(double **H, int n);
};