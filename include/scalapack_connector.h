#ifndef SCALAPACK_CONNECTOR_H
#define SCALAPACK_CONNECTOR_H

extern "C"
{
	void Cblacs_pinfo(int *myid, int *nprocs);
	void Cblacs_get(int icontxt, int what, int *val);
	void Cblacs_gridinit(int* context, char* order, int np_row, int np_col);
	void Cblacs_gridinfo(int icontxt, int* nprow, int *npcol, int *myprow, int *mypcol);
	
	void descinit_( 
		int *desc, const int *m, const int *n, const int *mb, const int *nb, 
		const int *irsrc, const int *icsrc, const int *ictxt,
		const int *lld, int *info);

	int numroc_( const int *n, const int *nb, const int *iproc, const int *srcproc, const int *nprocs );

	void pdgemv_(
		const char *transa,
		const int *M, const int *N,
		const double *alpha,
		const double *A, const int *IA, const int *JA, const int *DESCA,
		const double *B, const int *IB, const int *JB, const int *DESCB, const int *K, 
		const double *beta, double *C, const int *IC, const int *JC, const int *DESCC,const int *L);
	void pdsyevd_(
		char *jobz,char *uplo,
		int *n,double *a,int *ia,int *ja,
		int *desca,double *w,double *z,int *iz,int *jz,int *descz,
		double *work,int *lwork,int *iwork,int *liwork,int *info);
    void pdsyev_(
		char *jobz,char *uplo,
		int *n,double *a,int *ia,int *ja,
		int *desca,double *w,double *z,int *iz,int *jz,int *descz,
		double *work,int *lwork,int *info);
	void Cblacs_exit(int Contxt);
}
#endif