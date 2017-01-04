#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void currentShareY(Domain *D)
{



}

/*//
void MPI_TransferF_XplusFilter(Domain *D)
{
    int i,rank,numberData;
    int myrank, nTasks; 
    double *behindF;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=5;
    behindF=(double *)malloc(numberData*sizeof(double )); 
    
    //Transferring even ~ odd cores             
    behindF[0]=D->field[D->nx+2].E1;
    behindF[1]=D->field[D->nx+2].Pr;
    behindF[2]=D->field[D->nx+2].Pl;
    behindF[3]=D->field[D->nx+2].Sr;
    behindF[4]=D->field[D->nx+2].Sl;
        
    if(myrank%2==1)
    {
       MPI_Recv(behindF,numberData, MPI_DOUBLE, D->beforeCore, D->beforeCore, MPI_COMM_WORLD,&status);  
       D->field[0].E1=behindF[0];
       D->field[0].Pr=behindF[1];
       D->field[0].Pl=behindF[2];
       D->field[0].Sr=behindF[3];
       D->field[0].Sl=behindF[4];
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
       MPI_Send(behindF,numberData, MPI_DOUBLE, D->nextCore, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    behindF[0]=D->field[D->nx+2].E1;
    behindF[1]=D->field[D->nx+2].Pr;
    behindF[2]=D->field[D->nx+2].Pl;
    behindF[3]=D->field[D->nx+2].Sr;
    behindF[4]=D->field[D->nx+2].Sl;
        
    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(behindF,numberData, MPI_DOUBLE, D->beforeCore, D->beforeCore, MPI_COMM_WORLD,&status);  
       D->field[0].E1=behindF[0];
       D->field[0].Pr=behindF[1];
       D->field[0].Pl=behindF[2];
       D->field[0].Sr=behindF[3];
       D->field[0].Sl=behindF[4];
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
       MPI_Send(behindF,numberData, MPI_DOUBLE, D->nextCore, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(behindF);
}
*/


void MPI_TransferJ_DSX_Yplus(Domain *D)
{
    int i,j,k,numberData,start,end,nx,nySub,nzSub,ibegin,ibottom;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks, rank; 

    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    nx=D->nx;
    nySub=D->nySub;
    nzSub=D->nzSub;
    //first 3 is 3rd interpolation, 2nd 3 is Jx,Jy,Jz;
    numberData=3*3*(nx+5)*(nzSub+5);
    ibegin=0;
    ibottom=nx+5;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    rank=myrank%D->M;

    //Transferring even ~ odd cores 
    start=0; 
    for(j=0; j<3; j++)
      for(k=0; k<nzSub+5; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->YplusJ[start+i]=D->Jx[i][jend+j][k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->YplusJ[start+i]=D->Jy[i][jend+j][k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->YplusJ[start+i]=D->Jz[i][jend+j][k];
        start+=ibottom;
      }
        
    if(rank%2==1)
    {
       MPI_Recv(D->YplusJ,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=0; j<3; j++)
         for(k=0; k<nzSub+5; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             D->Jx[i][jstart+j][k]+=D->YplusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jy[i][jstart+j][k]+=D->YplusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jz[i][jstart+j][k]+=D->YplusJ[i+start];
           start+=ibottom;
         }
    }
    else if(rank%2==0 && rank!=D->M-1)
    {
      MPI_Send(D->YplusJ,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(j=0; j<3; j++)
      for(k=0; k<nzSub+5; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->YplusJ[start+i]=D->Jx[i][jend+j][k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->YplusJ[start+i]=D->Jy[i][jend+j][k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->YplusJ[start+i]=D->Jz[i][jend+j][k];
        start+=ibottom;
      }        
    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(D->YplusJ,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=0; j<3; j++)
         for(k=0; k<nzSub+5; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             D->Jx[i][jstart+j][k]+=D->YplusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jy[i][jstart+j][k]+=D->YplusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jz[i][jstart+j][k]+=D->YplusJ[i+start];
           start+=ibottom;
         }
    }
    else if(rank%2==1 && rank!=D->M-1)
    {
      MPI_Send(D->YplusJ,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    }     
    MPI_Barrier(MPI_COMM_WORLD);        
}


void MPI_TransferJ_DSX_Yminus(Domain *D)
{
    int i,j,k,numberData,start,end,nx,nySub,nzSub,ibegin,ibottom;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks, rank,share; 

    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    nx=D->nx;
    nySub=D->nySub;
    nzSub=D->nzSub;
    //first 2 is 3rd interpolation, 2nd 3 is Jx,Jy,Jz.
    numberData=2*3*(nx+5)*(nzSub+5);
    ibegin=0;
    ibottom=nx+5;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    rank=myrank%D->M;

    //Transferring even ~ odd cores 
    start=0; 
    for(j=1; j<3; j++)
      for(k=0; k<nzSub+5; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->YminusJ[start+i]=D->Jx[i][jstart-j][k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->YminusJ[start+i]=D->Jy[i][jstart-j][k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->YminusJ[start+i]=D->Jz[i][jstart-j][k];
        start+=ibottom;
      }
    if(rank%2==1 && rank!=D->M-1)
    {
       MPI_Recv(D->YminusJ,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=1; j<3; j++)
         for(k=0; k<nzSub+5; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             D->Jx[i][jend-j][k]+=D->YminusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jy[i][jend-j][k]+=D->YminusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jz[i][jend-j][k]+=D->YminusJ[i+start];
           start+=ibottom;
         }
    }
    else if(rank%2==0 && rank!=0)
    {
        
      MPI_Send(D->YminusJ,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(j=1; j<3; j++)
      for(k=0; k<nzSub+5; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->YminusJ[start+i]=D->Jx[i][jstart-j][k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->YminusJ[start+i]=D->Jy[i][jstart-j][k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->YminusJ[start+i]=D->Jz[i][jstart-j][k];
        start+=ibottom;
      }
        
    if(rank%2==0 && rank!=D->M-1)
    {
       MPI_Recv(D->YminusJ,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=1; j<3; j++)
         for(k=0; k<nzSub+5; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             D->Jx[i][jend-j][k]+=D->YminusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jy[i][jend-j][k]+=D->YminusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jz[i][jend-j][k]+=D->YminusJ[i+start];
           start+=ibottom;
         }
    }
    else if(rank%2==1)
    {
      MPI_Send(D->YminusJ,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_TransferF_DSX_YminusC(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int nz,int share)
{
    int i,j,k,numberData,start,end,ibegin,ibottom;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 

    MPI_Status status;         

    ibegin=0;
    ibottom=nx;   

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    //first 1 is 1 layer, 2nd 3 is 3 field variables, 3rd 1 is z.
    numberData=1*3*nx*nz; 
    rank=myrank%D->M;   

    //Transferring even ~ odd cores 
    start=0; 
    for(j=0; j<1; j++)	// 3 layers in y direction
      for(k=0; k<nz; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->minusYC[start+i]=f1[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusYC[start+i]=f2[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusYC[start+i]=f3[i][jstart+j][k];
        start+=nx;
      }

    if(rank%2==0 && rank!=D->M-1)
    {
      MPI_Recv(D->minusYC,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
      start=0; 
      for(j=0; j<1; j++)	// 3 layers in y direction
        for(k=0; k<nz; k++)
        {
          for(i=ibegin; i<ibottom; i++)
            f1[i][jend+j][k]=D->minusYC[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f2[i][jend+j][k]=D->minusYC[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f3[i][jend+j][k]=D->minusYC[start+i];
          start+=nx;
        }  
    }
    else if(rank%2==1)
       MPI_Send(D->minusYC,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(j=0; j<1; j++)	// 3 layers in y direction
      for(k=0; k<nz; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->minusYC[start+i]=f1[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusYC[start+i]=f2[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusYC[start+i]=f3[i][jstart+j][k];
        start+=nx;
      }
        
    if(rank%2==1 && rank!=D->M-1)
    {
      MPI_Recv(D->minusYC,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
      start=0; 
      for(j=0; j<1; j++)	// 1 layers in y direction
        for(k=0; k<nz; k++)
        {
          for(i=ibegin; i<ibottom; i++)
            f1[i][jend+j][k]=D->minusYC[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f2[i][jend+j][k]=D->minusYC[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f3[i][jend+j][k]=D->minusYC[start+i];
          start+=nx;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(D->minusYC,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

}

void MPI_TransferF_DSX_Yminus(Domain *D
         ,double ***f1,double ***f2,double ***f3
         ,double ***f4,double ***f5,double ***f6,int nx,int nz,int share)
{
    int i,j,k,numberData,start,end,ibegin,ibottom;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 

    MPI_Status status;         

    ibegin=0;
    ibottom=nx;   

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    //first 3 is 3 layer, 2nd 3 is 6 field variables
    numberData=share*6*nx*nz; 
    rank=myrank%D->M;   

    //Transferring even ~ odd cores 
    start=0; 
    for(j=0; j<share; j++)	// 3 layers in y direction
      for(k=0; k<nz; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f1[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f2[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f3[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f4[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f5[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f6[i][jstart+j][k];
        start+=nx;
      }

    if(rank%2==0 && rank!=D->M-1)
    {
      MPI_Recv(D->minusY,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
      start=0; 
      for(j=0; j<share; j++)	// 3 layers in y direction
        for(k=0; k<nz; k++)
        {
          for(i=ibegin; i<ibottom; i++)
            f1[i][jend+j][k]=D->minusY[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f2[i][jend+j][k]=D->minusY[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f3[i][jend+j][k]=D->minusY[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f4[i][jend+j][k]=D->minusY[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f5[i][jend+j][k]=D->minusY[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f6[i][jend+j][k]=D->minusY[start+i];
          start+=nx;
        }  
    }
    else if(rank%2==1)
       MPI_Send(D->minusY,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(j=0; j<share; j++)	// 3 layers in y direction
      for(k=0; k<nz; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f1[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f2[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f3[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f4[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f5[i][jstart+j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f6[i][jstart+j][k];
        start+=nx;
      }
        
    if(rank%2==1 && rank!=D->M-1)
    {
      MPI_Recv(D->minusY,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
      start=0; 
      for(j=0; j<share; j++)	// 1 layers in y direction
        for(k=0; k<nz; k++)
        {
          for(i=ibegin; i<ibottom; i++)
            f1[i][jend+j][k]=D->minusY[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f2[i][jend+j][k]=D->minusY[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f3[i][jend+j][k]=D->minusY[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f4[i][jend+j][k]=D->minusY[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f5[i][jend+j][k]=D->minusY[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f6[i][jend+j][k]=D->minusY[start+i];
          start+=nx;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(D->minusY,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

}

void MPI_TransferF_DSX_YplusC(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int nz,int share)
{
    int i,j,k,numberData,start,end,ibegin,ibottom;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 

    MPI_Status status;         
   
    ibegin=0;
    ibottom=nx;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    //first 1 is 1 layer, 2nd 3 is 3 field variables
    numberData=1*3*nx*nz; 
    rank=myrank%D->M;   

    //Transferring even ~ odd cores 
    start=0; 
    for(j=1; j<2; j++)
      for(k=0; k<nz; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->plusYC[start+i]=f1[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusYC[start+i]=f2[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusYC[start+i]=f3[i][jend-j][k];
        start+=nx;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(D->plusYC,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=1; j<2; j++)
         for(k=0; k<nz; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             f1[i][jstart-j][k]=D->plusYC[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f2[i][jstart-j][k]=D->plusYC[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f3[i][jstart-j][k]=D->plusYC[start+i];
           start+=nx;
         }
    }
    else if(rank%2==0 && rank!=D->M-1)
       MPI_Send(D->plusYC,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(j=1; j<2; j++)
      for(k=0; k<nz; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->plusYC[start+i]=f1[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusYC[start+i]=f2[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusYC[start+i]=f3[i][jend-j][k];
        start+=nx;
      }
        
    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(D->plusYC,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=1; j<2; j++)
         for(k=0; k<nz; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             f1[i][jstart-j][k]=D->plusYC[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f2[i][jstart-j][k]=D->plusYC[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f3[i][jstart-j][k]=D->plusYC[start+i];
           start+=nx;
         }
    }
    else if(rank%2==1 && rank!=D->M-1)
       MPI_Send(D->plusYC,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
}


void MPI_TransferF_DSX_Yplus(Domain *D
           ,double ***f1,double ***f2,double ***f3
           ,double ***f4,double ***f5,double ***f6,int nx,int nz,int share)
{
    int i,j,k,numberData,end,ibegin,ibottom;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank,start; 

    MPI_Status status;         
   
    ibegin=0;
    ibottom=nx;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    //first 2 is 2 layer, 2nd 6 is 6 field variables
    numberData=(share-1)*6*nx*nz; 
    rank=myrank%D->M;   

    //Transferring even ~ odd cores 
    start=0; 
    for(j=1; j<share; j++)
      for(k=0; k<nz; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f1[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f2[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f3[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f4[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f5[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f6[i][jend-j][k];
        start+=nx;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(D->plusY,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=1; j<share; j++)
         for(k=0; k<nz; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             f1[i][jstart-j][k]=D->plusY[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f2[i][jstart-j][k]=D->plusY[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f3[i][jstart-j][k]=D->plusY[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f4[i][jstart-j][k]=D->plusY[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f5[i][jstart-j][k]=D->plusY[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f6[i][jstart-j][k]=D->plusY[start+i];
           start+=nx;
         }
    }
    else if(rank%2==0 && rank!=D->M-1)
       MPI_Send(D->plusY,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(j=1; j<share; j++)
      for(k=0; k<nz; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f1[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f2[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f3[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f4[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f5[i][jend-j][k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f6[i][jend-j][k];
        start+=nx;
      }
        
    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(D->plusY,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=1; j<share; j++)
         for(k=0; k<nz; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             f1[i][jstart-j][k]=D->plusY[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f2[i][jstart-j][k]=D->plusY[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f3[i][jstart-j][k]=D->plusY[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f4[i][jstart-j][k]=D->plusY[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f5[i][jstart-j][k]=D->plusY[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f6[i][jstart-j][k]=D->plusY[start+i];
           start+=nx;
         }
    }
    else if(rank%2==1 && rank!=D->M-1)
       MPI_Send(D->plusY,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
}
