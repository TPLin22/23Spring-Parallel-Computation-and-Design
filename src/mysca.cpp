#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <lapacke.h>
#include <mpi.h>
#include "mysca.h"
#include "scalapack_connector.h"
using namespace std;

mysca::mysca(double ** H,int M,int myrank_mpi,int nprocs_mpi){
        ofstream ofs;
        int ictxt,nprow,npcol,myrow,mycol,mb,nb;
        int info,itemp;
        int ZERO = 0, ONE = 1;
        nprow = npcol = 2;
        mb = nb = 1;

        Cblacs_pinfo(&myrank_mpi, &nprocs_mpi);//tell scalapack size and rank
        Cblacs_get(-1, 0, &ictxt);//initialize ictxt
        Cblacs_gridinit(&ictxt, "Row", nprow, npcol);//tell sca nprow and npcol
        Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);//get my row and my col
   
        int mA = numroc_(&M,&nb,&myrow,&ZERO,&nprow);
        int nA = numroc_(&M,&nb,&mycol,&ZERO,&npcol);
        int mZ=mA,nZ=nA;
        int descA[10],descZ[10];

        descinit_(descA,&M,&M,&mb,&nb,&ZERO,&ZERO,&ictxt,&mA,&info);
        if (info != 0) {
            printf("Error in descinit_ for desca: %d\n", info);
            //return 1;
        }
        descinit_(descZ,&M,&M,&mb,&nb,&ZERO,&ZERO,&ictxt,&mZ,&info);
        if (info != 0) {
            printf("Error in descinit_ for desca: %d\n", info);
            //return 1;
        }
        A = new double[mA*nA];
        Z = new double[mZ*nZ];

        mA_array = new int[nprocs_mpi];
        nA_array = new int[nprocs_mpi];
        myrow_array = new int[nprocs_mpi];
        mycol_array = new int[nprocs_mpi];

        MPI_Allgather(&mA, 1, MPI_INT, mA_array, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&nA, 1, MPI_INT, nA_array, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&myrow, 1, MPI_INT, myrow_array, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&mycol, 1, MPI_INT, mycol_array, 1, MPI_INT, MPI_COMM_WORLD);

        hasm = new int[nprow],hasn = new int[npcol];
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
                //cout << "myrank=" << myrank_mpi << " ma_offset=" << mA_offset << " na_offset=" << nA_offset << " global i =" << global_i << " global j =" << global_j <<endl; 
                A[i*nA+j] = H[global_i][global_j];
            }
        }

        w = new double[M];
        work = new double[1];
        int lwork = -1;
        pdsyev_("V","U",&M,A,&ONE,&ONE,descA,w,Z,&ONE,&ONE,descZ,work,&lwork,&info);
        if (info != 0) {
            printf("Error in pdsyev_ for computing workspace size: %d\n", info);
            //return 1;
        }
        lwork = (int)work[0];
        delete[] work;
        work = new double[lwork];
        pdsyev_("V","U",&M,A,&ONE,&ONE,descA,w,Z,&ONE,&ONE,descZ,work,&lwork,&info);
        if (info != 0) {
            printf("Error in pdsyev_ for computing workspace size: %d\n", info);
            //return 1;
        }
        if (myrank_mpi == 0){
            ofs.open("../output/scaeigenvalues.log",ios::out);
            ofs << "using SCALAPACK:" << endl;
            for (int i = 0 ; i < M ; i ++)
            ofs << " w["  << i << "]=" << w[i] << endl;
            ofs.close();
        }

        GZ = new double[M*M];
        for (int i = 0 ; i < M*M ; i ++) GZ[i] = 0;
        for (int i = 0 ; i < mA ; i ++) {
            for (int j = 0 ; j < nA ; j ++) {
                int global_i = mA_offset + i;
                int global_j = nA_offset + j;
        //cout << "myrank=" << myrank_mpi << " ma_offset=" << mA_offset << " na_offset=" << nA_offset << " global i =" << global_i << " global j =" << global_j <<endl; 
        //A[i*nA+j] = H[global_i][global_j];
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

};
mysca::~mysca(){
        delete[] A;delete[] Z;delete[] work;delete[] w;
        delete[] mA_array;delete[] nA_array;delete[] myrow_array;delete[] mycol_array;
        delete[] hasm;delete[] hasn;
        delete[] GZ;
        A = NULL;Z = NULL;work = NULL;w = NULL;
        mA_array = NULL;nA_array = NULL;myrow_array = NULL;mycol_array = NULL;
        hasm = NULL; hasn = NULL;
        GZ = NULL;
}