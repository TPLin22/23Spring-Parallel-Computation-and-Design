#include <iostream>
#include <cstdio>
#include "Input_0.h"
#include "Input_p.h"
#include "Input_v.h"
#include "Input_d.h"
#include "mathwork.h"
#include <fstream>
#include <omp.h>
#include <iomanip>
#include <cmath>
#include <lapacke.h>
#include <mpi.h>
#include "scalapack_connector.h"
//#include "mysca.h"
using namespace std;

double** H;double* h;
int isSCALAPACK,M;

//calculate the position of every dV
void getlocation(double &x,double &y,double &z,double &dv,int i,int j,int k,int nx,int ny,int nz,double lx,double ly,double lz){
    x = y = z = 0.0;
    x = ((double)i/(double)nx)*lx;
    y = ((double)j/(double)ny)*ly;
    z = ((double)k/(double)nz)*lz;
    dv = (lx/(double)nx)*(ly/(double)ny)*(lz/(double)nz);
}

void Matirx_H_Create(int n){
    H = new double*[n];
    //#pragma omp parallel
    //{
        //#pragma omp for
        for (int i = 0 ; i < n ; i ++){
            H[i] = new double[n];
        }
        //#pragma omp for
        for (int i = 0 ; i < n ; i ++){
            for (int j = 0 ; j < n ; j ++)
                H[i][j] = 0;
        }
    //}
}

int main(int argc, char* argv[])
{
    int myrank_mpi, nprocs_mpi;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
    isSCALAPACK = 0;
    ofstream ofs;
    double start_time = MPI_Wtime(),Pre_time,Now_time;
if (myrank_mpi == 0)
{
//read maininput
    Input_0 ourInput = Input_0(argv[1]);
//ourInput store path of files,the length of x,y,z(lx,ly,lz)

//read the position of point
    Input_p ourPos = Input_p(ourInput.points_path.data());
    M = ourPos.n;
//ourPos store the num(n) and location(pos[i][1,2,3]) of the points,

//read the venergy of points
    Input_v ourV = Input_v(ourInput.v_path.data());
//ourV store the num of grid(nx,ny,nz) and venergy of each grid(v[x][y][z])
    
//read the distribution function
    Input_d ourD = Input_d(ourInput.distribution_path.data());
//ourD store the cutoff radius ,dr and values of f(r)

//initialize something of f(r),used for interpolations
    mathwork ourM = mathwork(ourD.f,ourD.mesh,ourD.dr);

    int nthreads = 16;
    omp_set_num_threads(nthreads);

//create matrix H
    Matirx_H_Create(ourPos.n);
    Now_time = MPI_Wtime();
    cout.fill('-');
    cout << setw(20) << "OPERATION" << setw(20) << "TIME(sec)" <<endl;
    cout.fill(' ');
    cout << setw(20) << "INITIALIZE" << setw(20) << Now_time-start_time << endl;
    Pre_time = MPI_Wtime();

//if size of point is big,do some optimization
    
//calculate integral 
    //#pragma omp parallel
    //{
        /*#pragma omp for
        for (int i = 0 ; i < ourV.nx ; i ++){
            int numpos = 0;
            int Vpos[50],pos1,pos2;
            double x,y,z,dv,r,vans;
            for (int j = 0 ; j < ourV.ny ; j ++)
                for (int k = 0 ; k < ourV.nz ; k ++){
                    
                    getlocation(x,y,z,dv,i,j,k,ourV.nx,ourV.ny,ourV.nz,ourInput.lx,ourInput.ly,ourInput.lz);
                    numpos = 0;
                    //search which points are close enough to this grid(r < cutoff),store them in Vpos
                    for (int l = 0 ; l < M ; l ++){
                    r = sqrt(pow(x-ourPos.pos[l][1],2)+pow(y-ourPos.pos[l][2],2)+pow(z-ourPos.pos[l][3],2));
                    if (r <= ourD.cutoff) Vpos[numpos++] = l;
                    }
                    if (numpos == 0) continue;
                    //ergodic these points,calculate the integral
                    for (int i2 = 0 ; i2 < numpos ; i2 ++)
                        for (int j2 = i2 ; j2 < numpos ; j2 ++){
                            pos1 = Vpos[i2],pos2 = Vpos[j2];
                            vans = ourM.calculat(x,y,z,ourPos.pos[pos1],ourPos.pos[pos2],ourV.venergy[i][j][k]);
                            #pragma omp atomic
                            H[pos1][pos2] += (vans*dv);
                            if (pos1 != pos2)
                            {
                            #pragma omp atomic
                            H[pos2][pos1] += (vans*dv);
                            }
                        }
                }
        }*/
        double dwx = ourInput.lx / (double)ourV.nx;
        double dwy = ourInput.ly / (double)ourV.ny;
        double dwz = ourInput.lz / (double)ourV.nz;
        for (int i = 0 ; i < M ; i ++){
            for (int j = i ; j < M ; j ++){
                double r = sqrt(pow(ourPos.pos[j][1]-ourPos.pos[i][1],2)+pow(ourPos.pos[j][2]-ourPos.pos[i][2],2)+pow(ourPos.pos[j][3]-ourPos.pos[i][3],2));
                if (r > ourD.cutoff) continue;
                //cout << "i = " << i << " j = " << j << endl;
                int ibegin = max((int)((ourPos.pos[i][1]-ourD.cutoff)/dwx),(int)((ourPos.pos[j][1]-ourD.cutoff)/dwx));ibegin = max(ibegin,0);
                int iend = min((int)((ourPos.pos[i][1]+ourD.cutoff)/dwx),(int)((ourPos.pos[j][1]+ourD.cutoff)/dwx));iend = min(iend,ourV.nx);
                int jbegin = max((int)((ourPos.pos[i][2]-ourD.cutoff)/dwy),(int)((ourPos.pos[j][2]-ourD.cutoff)/dwy));jbegin = max(jbegin,0);
                int jend = min((int)((ourPos.pos[i][2]+ourD.cutoff)/dwy),(int)((ourPos.pos[j][2]+ourD.cutoff)/dwy));jend = min(jend,ourV.ny);
                int kbegin = max((int)((ourPos.pos[i][3]-ourD.cutoff)/dwz),(int)((ourPos.pos[j][3]-ourD.cutoff)/dwz));kbegin = max(kbegin,0);
                int kend = min((int)((ourPos.pos[i][3]+ourD.cutoff)/dwz),(int)((ourPos.pos[j][3]+ourD.cutoff)/dwz));kend = min(kend,ourV.nz);
                /*cout << "1x = " << ourPos.pos[i][1] << " 2x = " << ourPos.pos[j][1] << endl;
                cout << "1left = " << ourPos.pos[i][1]-ourD.cutoff << " 2left =  " << ourPos.pos[j][1]-ourD.cutoff << endl;
                cout << "ibegin = " << ibegin << " iend = " << iend << endl;
                cout << "jbegin = " << jbegin << " jend = " << jend << endl;
                cout << "kbegin = " << kbegin << " kend = " << kend << endl;*/
                #pragma omp parallel for
                for (int i1 = ibegin ; i1 < iend ; i1 ++){
                    for (int j1 = jbegin ; j1 < jend ; j1 ++)
                        for (int k1 = kbegin ; k1 < kend ; k1 ++){
                            double x,y,z,dv,vans;
                            getlocation(x,y,z,dv,i1,j1,k1,ourV.nx,ourV.ny,ourV.nz,ourInput.lx,ourInput.ly,ourInput.lz);
                            if (i1 < 0 || j1 < 0 || k1 < 0 || i1 > ourV.nx || j1 > ourV.ny || k1 > ourV.nz) cout << "error!" << endl;
                            vans = ourM.calculat(x,y,z,ourPos.pos[i],ourPos.pos[j],ourV.venergy[i1][j1][k1]);
                            #pragma omp atomic
                            H[i][j] += (vans*dv);
                            if (i != j)
                            {
                            #pragma omp atomic
                            H[j][i] += (vans*dv);
                            }
                        }
                }
            }
        }
    //}

//counting time using for calculating integral
    Now_time = MPI_Wtime();
    cout << setw(20) << "INTEGRAL" << setw(20) << Now_time-Pre_time << endl;
    Pre_time = MPI_Wtime();

//outprint the matrix H to out.log
    ofs.open("../output/out.log",ios::out);
    for (int i = 0 ; i < ourPos.n ; i ++){
        for (int j = 0 ; j <ourPos.n ; j ++)
            ofs << setw(10) << setprecision(8) << H[i][j] << " ";ofs << endl;
    }
    ofs.close();

//use lapack to diagonize
    if (ourInput.diago_lib == "lapack"){
        ourM.uselapack(H,ourPos.n);
        for (int i = 0 ; i < ourPos.n ; i ++){
            delete[] H[i];
            H[i] = NULL;
        }
        delete[] H;//release matrix H
        H = NULL;
        Now_time = MPI_Wtime();
        cout << setw(20) << "LAPACK" << setw(20) << Now_time-Pre_time << endl;
        cout << setw(20) << "TOTAL" << setw(20) << Now_time-start_time << endl;
    }
    else{//mark isSCALAPACK to let other process know to call scalapack
        isSCALAPACK = 1;
        h = new double[M*M];//use one dimension array to more easily send messages using MPI
        for (int i = 0 ; i < M ; i ++)
            for (int j = 0 ; j < M ; j ++)
                h[i*M+j] = H[i][j];
    }
}
//Bcast whether to call scalapack or not
    MPI_Bcast(&isSCALAPACK,1,MPI_INT,0,MPI_COMM_WORLD);
//using scalapack to diagonize
    if (isSCALAPACK == 1){
        //Bcast M(dimension of global Matrix H and the elements of H)
        MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (myrank_mpi != 0){//rank 0 has already created the matrix
            h = new double[M*M];
        }
        MPI_Bcast(h, M * M, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        MPI_Barrier(MPI_COMM_WORLD);

        if (myrank_mpi != 0){
            H = new double*[M];
            for (int i = 0 ; i <M ; i ++)
                H[i] = new double[M];
            for (int i = 0 ; i < M ; i ++)
                for (int j = 0 ; j <M ; j ++)
                    H[i][j] = h[i*M+j];
        }
        delete[] h;
        h = NULL;

        int ictxt,nprow,npcol,myrow,mycol,mb,nb;
        int info,itemp;
        int ZERO = 0, ONE = 1;
        nprow = npcol = 2;
        mb = nb = 1;

        Cblacs_pinfo(&myrank_mpi, &nprocs_mpi);//tell scalapack size and rank
        Cblacs_get(-1, 0, &ictxt);//initialize ictxt
        Cblacs_gridinit(&ictxt, "Row", nprow, npcol);//tell scalapack nprow and npcol
        Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);//get my row and my col
   
        int mA = numroc_(&M,&nb,&myrow,&ZERO,&nprow);//get row width of this grid
        int nA = numroc_(&M,&nb,&mycol,&ZERO,&npcol);//get col width of this grid
        int mZ=mA,nZ=nA;
        int descA[10],descZ[10];

        descinit_(descA,&M,&M,&mb,&nb,&ZERO,&ZERO,&ictxt,&mA,&info);//initialize descA which store the messages of A
        if (info != 0) {
            printf("Error in descinit_ for desca: %d\n", info);
            return 1;
        }
        descinit_(descZ,&M,&M,&mb,&nb,&ZERO,&ZERO,&ictxt,&mZ,&info);//initialize descZ which store the messages of Z
        if (info != 0) {
            printf("Error in descinit_ for desca: %d\n", info);
            return 1;
        }
        double *A = new double[mA*nA];
        double *Z = new double[mZ*nZ];

        int* mA_array = new int[nprocs_mpi];
        int* nA_array = new int[nprocs_mpi];
        int* myrow_array = new int[nprocs_mpi];
        int* mycol_array = new int[nprocs_mpi];

//All process get the mA of nA of every process,which is used to calculate the global location of local matrix
        MPI_Allgather(&mA, 1, MPI_INT, mA_array, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&nA, 1, MPI_INT, nA_array, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&myrow, 1, MPI_INT, myrow_array, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&mycol, 1, MPI_INT, mycol_array, 1, MPI_INT, MPI_COMM_WORLD);

        int *hasm = new int[nprow],*hasn = new int[npcol];
        int mA_offset = 0;
        for (int i = 0 ; i < nprow ; i ++) hasm[i] = 0;
        for (int j = 0 ; j < npcol ; j ++) hasn[j] = 0;
        for (int p = 0; p < myrank_mpi; ++p)
            if (myrow_array[p] < myrow && hasm[myrow_array[p]] == 0){
                mA_offset += mA_array[p];
                hasm[myrow_array[p]] = 1;
            }
        int nA_offset = 0;
        for (int p = 0; p < myrank_mpi; ++p)
            if (mycol_array[p] < mycol && hasn[mycol_array[p]] == 0){
                nA_offset += nA_array[p];
                hasn[mycol_array[p]] = 1;
            }

        for (int i = 0 ; i < mA ; i ++) {
            for (int j = 0 ; j < nA ; j ++) {
                int global_i = mA_offset + i;
                int global_j = nA_offset + j;
                A[i*nA+j] = H[global_i][global_j];
            }
        }
        double *w = new double[M];
        double *work = new double[1];
        int lwork = -1;//set lwork = -1 to query work space
        pdsyev_("V","U",&M,A,&ONE,&ONE,descA,w,Z,&ONE,&ONE,descZ,work,&lwork,&info);
        if (info != 0) {
            printf("Error in pdsyev_ for computing workspace size: %d\n", info);
            return 1;
        }
        lwork = (int)work[0];
        delete[] work;
        work = new double[lwork];
        pdsyev_("V","U",&M,A,&ONE,&ONE,descA,w,Z,&ONE,&ONE,descZ,work,&lwork,&info);//call pdsyev to diagonize it 
        if (info != 0) {
            printf("Error in pdsyev_ for computing workspace size: %d\n", info);
            return 1;
        }
        if (myrank_mpi == 0){
            ofs.open("../output/scaeigenvalues.log",ios::out);
            ofs << "using SCALAPACK:" << endl;
            for (int i = 0 ; i < M ; i ++)
            ofs << " w["  << i << "]=" << w[i] << endl;
            ofs.close();
        }
//since Z is local matrix,every process should send Z back to process 0 to make the global Matrix GZ
        double *GZ = new double[M*M];
        for (int i = 0 ; i < M*M ; i ++) GZ[i] = 0;
        for (int i = 0 ; i < mA ; i ++) {
            for (int j = 0 ; j < nA ; j ++) {
                int global_i = mA_offset + i;
                int global_j = nA_offset + j;
                GZ[global_i*M+global_j] = Z[i*nA+j];
            }
        }
        if (myrank_mpi == 0) {
            double* recvbuf = new double[M*M];
            MPI_Reduce(GZ, recvbuf, M*M, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            for (int i = 0 ; i < M*M ; i ++) GZ[i] = recvbuf[i];
            delete[] recvbuf;
        } else {
            MPI_Reduce(GZ, NULL, M*M, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        if (myrank_mpi == 0) {
            ofs.open("../output/scaeigenvectors.log",ios::out);
            ofs << "using SCALAPACK:" << endl;
            for (int i = 0; i < M; i++) {
                ofs << "v[" << i << "] = ";
                for (int j = 0; j < M; j++) {
                    ofs << GZ[i * M + j] << " ";
                }
                 ofs << endl;
            }
            ofs.close();
        }
        delete[] A;delete[] Z;delete[] work;delete[] w;
        delete[] mA_array;delete[] nA_array;delete[] myrow_array;delete[] mycol_array;
        delete[] hasm;delete[] hasn;
        delete[] GZ;
        A = NULL;Z = NULL;work = NULL;w = NULL;
        mA_array = NULL;nA_array = NULL;myrow_array = NULL;mycol_array = NULL;
        hasm = NULL; hasn = NULL;
        GZ = NULL;
        for (int i = 0 ; i < M ; i ++){
            delete[] H[i];
            H[i] = NULL;
        }
        delete[] H;
        H = NULL;
        if (myrank_mpi == 0)
        {Now_time = MPI_Wtime();
        cout << setw(20) << "SCALAPACK" << setw(20) << Now_time-Pre_time << endl;
        cout << setw(20) << "TOTAL" << setw(20) << Now_time-start_time << endl;}
    }
    
    MPI_Finalize();
    return 0;
}

//
//build the spatical net,a n*n*n array
//it is ok to store f[r] with mesh double
//ergodic ngrid points,count 1 * v * 2 (O(ngrid))
//ergodic ngrid points, ergodic given points (50*50) (O(ngrid*2500)) 
// if ngrid = around 500^3, it needs 312,500,000,000 (1 hour)
// if ngrid = around 50^3, it needs 312,500,000 (3 seconds)

//first ergodic ngrid, delete the points which is not include in any ball?
//store the points in a new array, it can cut down several useless points
//even it can store what balls a point is in, it will only contribute to needed matrix elements
//may be it can cut down from 2500 to 1

