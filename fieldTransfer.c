#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "mesh.h"
#include "constants.h"

void MPI_TransferF_DSX_YminusC(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int ny,int share);
void MPI_TransferF_DSX_YplusC(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int ny,int share);
void MPI_TransferF_DSX_Yminus(Domain *D,
       double ***f1,double ***f2,double ***f3,
       double ***f4,double ***f5,double ***f6,int nx,int ny,int share);
void MPI_TransferF_DSX_Yplus(Domain *D,
       double ***f1,double ***f2,double ***f3,
       double ***f4,double ***f5,double ***f6,int nx,int ny,int share);
void MPI_TransferF_DSX_ZminusC(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int ny,int share);
void MPI_TransferF_DSX_ZplusC(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int ny,int share);
void MPI_TransferF_DSX_Zminus(Domain *D,
       double ***f1,double ***f2,double ***f3,
       double ***f4,double ***f5,double ***f6,int nx,int ny,int share);
void MPI_TransferF_DSX_Zplus(Domain *D,
       double ***f1,double ***f2,double ***f3,
       double ***f4,double ***f5,double ***f6,int nx,int ny,int share);

void fieldTransferC(Domain *D)
{
  int myrank, nTasks,rank,rankM,rankN;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch((D->fieldType-1)*3+D->dimension) {
  //1D field
  case (Split-1)*3+1:
    break;
  //split 2D
  case (Split-1)*3+2:
    //D->nx+5, 1 is number z, 3 is share.
    MPI_TransferF_DSX_YminusC(D,D->SrC,D->SlC,D->ExC,D->nx+5,1,3);
    MPI_TransferF_DSX_YplusC(D,D->PrC,D->PlC,D->BxC,D->nx+5,1,3);
    MPI_Barrier(MPI_COMM_WORLD);
    break;

  //3D
  case (Split-1)*3+3:
    //3 is share. But it works on '..C'.
    if(D->M>1)
    {
      MPI_TransferF_DSX_YminusC(D,D->SrC,D->SlC,D->ExC,D->nx+5,D->nzSub+5,3);
      MPI_TransferF_DSX_YplusC(D,D->PrC,D->PlC,D->BxC,D->nx+5,D->nzSub+5,3);
    }
    else	;
    if(D->N>1)
    {
      MPI_TransferF_DSX_ZminusC(D,D->PrC,D->PlC,D->ExC,D->nx+5,D->nySub+5,3);
      MPI_TransferF_DSX_ZplusC(D,D->SrC,D->SlC,D->BxC,D->nx+5,D->nySub+5,3);
    }
    else	;
    MPI_Barrier(MPI_COMM_WORLD);
    break;
  default:
    printf("what fieldType? and what dimension?\n");
  }
}

void fieldTransfer(Domain *D)
{
  int myrank, nTasks,rank,rankM,rankN;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch((D->fieldType-1)*3+D->dimension) {
  //1D field
  case (Split-1)*3+1:
    break;
  //split 2D
  case (Split-1)*3+2:
    MPI_TransferF_DSX_Yminus(D,D->Sr,D->Sl,D->Ex,D->Pr,D->Pl,D->Bx,D->nx+5,1,3);
    MPI_TransferF_DSX_Yplus(D,D->Pr,D->Pl,D->Bx,D->Sr,D->Sl,D->Ex,D->nx+5,1,3);
    MPI_Barrier(MPI_COMM_WORLD);
    break;

  //3D
  case (Split-1)*3+3:
    if(D->M>1)
    {
      MPI_TransferF_DSX_Yminus(D,D->Sr,D->Sl,D->Ex,D->Pr,D->Pl,D->Bx,D->nx+5,D->nzSub+5,3);
      MPI_TransferF_DSX_Yplus(D,D->Pr,D->Pl,D->Bx,D->Sr,D->Sl,D->Ex,D->nx+5,D->nzSub+5,3);
    }
    else	;
    if(D->N>1)
    {
      MPI_TransferF_DSX_Zminus(D,D->Sr,D->Sl,D->Ex,D->Pr,D->Pl,D->Bx,D->nx+5,D->nySub+5,3);
      MPI_TransferF_DSX_Zplus(D,D->Pr,D->Pl,D->Bx,D->Sr,D->Sl,D->Ex,D->nx+5,D->nySub+5,3);
    }
    else	;
    MPI_Barrier(MPI_COMM_WORLD);
    break;
  default:
    printf("what fieldType? and what dimension?\n");
  }
}

