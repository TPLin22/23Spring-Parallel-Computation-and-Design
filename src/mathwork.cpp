#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <lapacke.h>
#include "mathwork.h"
using namespace std;

mathwork::mathwork(double *f, int size_of_f, double dr){
//f store given distribution function,bellow array are used for preprocess
    size = size_of_f;
    xdata = new double[size];
    ydata = new double[size];
    h = new double[size];
    aver = new double[size];
    lan = new double[size];
    d = new double[size];
    nu = new double[size];
    a = new double*[size];
    for (int i = 0 ; i < size ; i ++)
        a[i] = new double[size];
    beita = new double[size];
    yy = new double[size];
    xx = new double[size];
    xdata[0] = 0;
    for (int i = 1 ; i <size ; i ++){
        xdata[i] = xdata[i-1] + dr;
        ydata[i] = f[i];
    }
    ydata[0] = f[0];


    for (int i=0;i<size-1;i++){
   	    h[i]=xdata[i+1]-xdata[i];
    } 
    for (int i=0;i<size-1;i++){
	    aver[i]=(ydata[i+1]-ydata[i])/h[i];		
	}
        lan[0]=1;
    	d[0]=(6/h[0])*(aver[0]-0);
    	nu[size-2]=1;
    	d[size-1]=(6/h[size-2])*(0-aver[size-2]);

    for (int i=0;i<size-2;i++){
		nu[i]=h[i]/(h[i]+h[i+1]);
	}
	for (int i=1;i<size-1;i++){
		lan[i]=h[i]/(h[i-1]+h[i]);
	}
	for (int i=1;i<size-1;i++){
		d[i]=6*(aver[i]-aver[i-1])/(h[i-1]+h[i]);
	}
    for (int i=0;i<size;i++){
	 	a[i][i]=2;
	 }
	 for (int i=0;i<size-1;i++){
	 	a[i+1][i]=nu[i];
	 }
	 for (int i=0;i<size-1;i++){
	 	a[i][i+1]=lan[i];
	 }	
     beita[0]=lan[0]/2;
     for (int i=1;i<size-1;i++){
    	beita[i]=lan[i]/(2-nu[i-1]*beita[i-1]);
	}
    yy[0]=d[0]/2;
    for (int i=1;i<size;i++){
    yy[i]=(d[i]-nu[i-1]*yy[i-1])/(2-nu[i-1]*beita[i-1]);	
	}
    xx[size-1]=yy[size-1];
	for (int i=size-2;i>=0;i--){
		xx[i]=yy[i]-beita[i]*xx[i+1];
	}
}

mathwork::~mathwork(){
    delete[] xdata;
    delete[] ydata;
    delete[] h;
    delete[] aver;
    delete[] lan;
    delete[] d;
    delete[] nu;
    delete[] beita;
    delete[] yy;
    delete[] xx;
    for (int i = 0 ; i < size ; i ++){
        delete a[i];
        a[i] = NULL;
    }
    a = NULL;
    xdata = NULL;
    ydata = NULL;
    h = NULL;
    aver = NULL;
    lan = NULL;
    d = NULL;
    nu = NULL;
    beita = NULL;
    yy = NULL;
    xx = NULL;
}

void mathwork::uselapack(double **H, int n){
    ofstream ofs;
        double *A;double *w;
        A = new double[n*n];
        w = new double[n];
        for (int i = 0 ; i < n ; i ++)
            for (int j = 0 ; j < n ; j ++)
                A[i*n+j] = H[i][j];
        int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',n,A,n,w);//use dsyev
        if (info != 0){
            cout << "FAILED TO DIAGONALIZE THIS MATRIX!" << endl;
        }
        else{
            ofs.open("../output/lapackeigenvalues.log",ios::out);
            //freopen("../output/lapackeigenvalues.log","w",stdout);
            ofs << "using LAPACK:" << endl;
            ofs << "Eigenvalues: ";//eigenvalues
            for (int i = 0 ; i < n ; i ++)
                ofs <<"w[" << i << "] = " << w[i] << endl;
            ofs << endl;
            ofs.close();
            ofs.open("../output/lapackeigenvectors.log",ios::out);
            //freopen("../output/lapackeigenvectors.log","w",stdout);
            ofs << "using LAPACK:" << endl;
            ofs << "Eigenvectors: " << endl;//eigenvectors
            for (int i = 0 ; i < n ; i ++){
                ofs << "v[" << i << "] = ";
                for (int j = 0 ; j < n ; j ++)
                    ofs << A[i*n+j] << " ";
                ofs << endl;
            }
            ofs.close();
        }
        delete[] A;
        A = NULL;
        delete[] w;
        w = NULL;
}

double mathwork::fl(double x){
//binary search the location to interpolate
    if (x > xdata[size-1]) return 0;
    int klo = 0;
    int khi = size - 1;
    while (khi - klo > 1){
        int k = (khi + klo) / 2;
        if (xdata[k] > x){
            khi = k;
        }  
        else{
            if (xdata[k] < x) klo = k;
            else return ydata[k];
        }
    }    
    int as = klo;

//interpolation
    if (as != size-1){
        double numble;
	    double sum1;
	    double sum2;
	    double sum3;
	    double sum4;
	    sum1=xx[as]*(xdata[as+1]-x)*(xdata[as+1]-x)*(xdata[as+1]-x)/6/h[as];
	    sum2=xx[as+1]*(x-xdata[as])*(x-xdata[as])*(x-xdata[as])/6/h[as];
	    sum3=(ydata[as]-xx[as]*h[as]*h[as]/6)*(xdata[as+1]-x)/h[as];
	    sum4=(ydata[as+1]-xx[as+1]*h[as]*h[as]/6)*(x-xdata[as])/h[as];
	    numble=sum1+sum2+sum3+sum4;
        return numble;
    }
    else return ydata[size-1];
}

double mathwork::calculat(double x,double y,double z,double x1[],double x2[],double v){
//calculate the integral
    double r1 = sqrt(pow(x-x1[1],2)+pow(y-x1[2],2)+pow(z-x1[3],2));
    double r2 = sqrt(pow(x-x2[1],2)+pow(y-x2[2],2)+pow(z-x2[3],2));
    double f1 = fl(r1);
    double f2 = fl(r2);
    return f1*f2*v;
}