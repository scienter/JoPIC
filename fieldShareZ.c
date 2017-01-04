#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>
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


void MPI_TransferJ_Zplus(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int ny)
{
    int i,k,j,numberData,start,end,ibegin,ibottom;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks, rank; 

    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    //first 3 is 3rd interpolation, 2nd 3 is Jx,Jy,Jz
    numberData=3*3*(nx)*(ny);
    ibegin=0;
    ibottom=nx;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    rank=(int)(myrank/D->M);

    //Transferring even ~ odd cores 
    start=0;
    for(k=0; k<3; k++)
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->ZplusJ[i+start]=D->Jx[i][j][kend+k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->ZplusJ[i+start]=D->Jy[i][j][kend+k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->ZplusJ[i+start]=D->Jz[i][j][kend+k];
        start+=ibottom;
      }
    if(rank%2==1)
    {
       MPI_Recv(D->ZplusJ,numberData, MPI_DOUBLE, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);  
       start=0;
       for(k=0; k<3; k++)
         for(j=0; j<ny; j++)
         {
           for(i=ibegin; i<ibottom; i++)
             D->Jx[i][j][kstart+k]+=D->ZplusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jy[i][j][kstart+k]+=D->ZplusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jz[i][j][kstart+k]+=D->ZplusJ[i+start];
           start+=ibottom;
         }
    }
    else if(rank%2==0 && rank!=D->N-1)
    {
      MPI_Send(D->ZplusJ,numberData, MPI_DOUBLE, myrank+D->M, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(k=0; k<3; k++)
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->ZplusJ[i+start]=D->Jx[i][j][kend+k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->ZplusJ[i+start]=D->Jy[i][j][kend+k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->ZplusJ[i+start]=D->Jz[i][j][kend+k];
        start+=ibottom;
      }
    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(D->ZplusJ,numberData, MPI_DOUBLE, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);  
       start=0;
       for(k=0; k<3; k++)
         for(j=0; j<ny; j++)
         {
           for(i=ibegin; i<ibottom; i++)
             D->Jx[i][j][kstart+k]+=D->ZplusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jy[i][j][kstart+k]+=D->ZplusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jz[i][j][kstart+k]+=D->ZplusJ[i+start];
           start+=ibottom;
         }
    }
    else if(rank%2==1 && rank!=D->N-1)
    {
      MPI_Send(D->ZplusJ,numberData, MPI_DOUBLE, myrank+D->M, myrank, MPI_COMM_WORLD);             
    }     
    MPI_Barrier(MPI_COMM_WORLD);        
}

void MPI_TransferJ_Zminus(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int ny)
{
    int i,j,k,numberData,start,end,ibegin,ibottom;
    int istart,iend,jstart,jend,kend,kstart;
    int myrank, nTasks,rank; 

    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    //first 2 is 3rd interpolation, 2nd 3 is Jx, Jy,Jz
    numberData=2*3*(nx)*(ny);
    ibegin=0;
    ibottom=nx;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    rank=(int)(myrank/D->M);

    //Transferring even ~ odd cores 
    start=0;
    for(k=1; k<3; k++)
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->ZminusJ[i+start]=D->Jx[i][j][kstart-k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->ZminusJ[i+start]=D->Jy[i][j][kstart-k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->ZminusJ[i+start]=D->Jz[i][j][kstart-k];
        start+=ibottom;
      }
    if(rank%2==1 && rank!=D->N-1)
    {
       MPI_Recv(D->ZminusJ,numberData, MPI_DOUBLE, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);  
       start=0;
       for(k=1; k<3; k++)
         for(j=0; j<ny; j++)
         {
           for(i=ibegin; i<ibottom; i++)
             D->Jx[i][j][kend-k]+=D->ZminusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jy[i][j][kend-k]+=D->ZminusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jz[i][j][kend-k]+=D->ZminusJ[i+start];
           start+=ibottom;
         }
    }
    else if(rank%2==0 && rank!=0)
    {
      MPI_Send(D->ZminusJ,numberData, MPI_DOUBLE, myrank-D->M, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(k=1; k<3; k++)
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->ZminusJ[i+start]=D->Jx[i][j][kstart-k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->ZminusJ[i+start]=D->Jy[i][j][kstart-k];
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->ZminusJ[i+start]=D->Jz[i][j][kstart-k];
        start+=ibottom;
      }
    if(rank%2==0 && rank!=D->N-1)
    {
       MPI_Recv(D->ZminusJ,numberData, MPI_DOUBLE, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);  
       start=0;
       for(k=1; k<3; k++)
         for(j=0; j<ny; j++)
         {
           for(i=ibegin; i<ibottom; i++)
             D->Jx[i][j][kend-k]+=D->ZminusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jy[i][j][kend-k]+=D->ZminusJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jz[i][j][kend-k]+=D->ZminusJ[i+start];
           start+=ibottom;
         }
    }
    else if(rank%2==1)
    {
      MPI_Send(D->ZminusJ,numberData, MPI_DOUBLE, myrank-D->M, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_TransferF_DSX_ZminusC(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int ny,int share)
{
    int i,j,k,numberData,start,end,nySub,nzSub,ibegin,ibottom;
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
    //first 1 is 1 layer, second 3 is 3 field variables
//    numberData=1*3*(nx+5)*(nySub+5); 
    rank=(int)(myrank/D->M);   

    //Transferring even ~ odd cores 
    start=0; 
    for(k=0; k<1; k++)	// 3 layers in z direction
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->minusZC[start+i]=f1[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZC[start+i]=f2[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZC[start+i]=f3[i][j][kstart+k];
        start+=nx;
      }

    if(rank%2==0 && rank!=D->N-1)
    {
      MPI_Recv(D->minusZC,D->numMinusZC, MPI_DOUBLE, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);  
      start=0; 
      for(k=0; k<1; k++)	// 3 layers in z direction
        for(j=0; j<ny; j++)
        {
          for(i=ibegin; i<ibottom; i++)
            f1[i][j][kend+k]=D->minusZC[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f2[i][j][kend+k]=D->minusZC[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f3[i][j][kend+k]=D->minusZC[start+i];
          start+=nx;
        }  
    }
    else if(rank%2==1)
       MPI_Send(D->minusZC,D->numMinusZC, MPI_DOUBLE, myrank-D->M, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(k=0; k<1; k++)	// 3 layers in z direction
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->minusZC[start+i]=f1[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZC[start+i]=f2[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZC[start+i]=f3[i][j][kstart+k];
        start+=nx;
      }
        
    if(rank%2==1 && rank!=D->N-1)
    {
      MPI_Recv(D->minusZC,D->numMinusZC, MPI_DOUBLE, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);  
      start=0; 
      for(k=0; k<1; k++)	// 3 layers in z direction
        for(j=0; j<ny; j++)
        {
          for(i=ibegin; i<ibottom; i++)
            f1[i][j][kend+k]=D->minusZC[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f2[i][j][kend+k]=D->minusZC[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f3[i][j][kend+k]=D->minusZC[start+i];
          start+=nx;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(D->minusZC,D->numMinusZC, MPI_DOUBLE, myrank-D->M, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

}

void MPI_TransferF_DSX_ZplusC(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int ny,int share)
{
    int i,j,k,numberData,start,end,nySub,nzSub,ibegin,ibottom;
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
//    numberData=1*3*(nx+5)*(nySub+5); 
    rank=(int)(myrank/D->M);   

    //Transferring even ~ odd cores 
    start=0; 
    for(k=1; k<2; k++)
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->plusZC[start+i]=f1[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZC[start+i]=f2[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZC[start+i]=f3[i][j][kend-k];
        start+=nx;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(D->plusZC,D->numPlusZC, MPI_DOUBLE, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);  
       start=0;
       for(k=1; k<2; k++)
         for(j=0; j<ny; j++)
         {
           for(i=ibegin; i<ibottom; i++)
             f1[i][j][kstart-k]=D->plusZC[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f2[i][j][kstart-k]=D->plusZC[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f3[i][j][kstart-k]=D->plusZC[start+i];
           start+=nx;
         }
    }
    else if(rank%2==0 && rank!=D->N-1)
       MPI_Send(D->plusZC,D->numPlusZC, MPI_DOUBLE, myrank+D->M, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(k=1; k<2; k++)
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->plusZC[start+i]=f1[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZC[start+i]=f2[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZC[start+i]=f3[i][j][kend-k];
        start+=nx;
      }
        
    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(D->plusZC,D->numPlusZC, MPI_DOUBLE, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);  
       start=0;
       for(k=1; k<2; k++)
         for(j=0; j<ny; j++)
         {
           for(i=ibegin; i<ibottom; i++)
             f1[i][j][kstart-k]=D->plusZC[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f2[i][j][kstart-k]=D->plusZC[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f3[i][j][kstart-k]=D->plusZC[start+i];
           start+=nx;
         }
    }
    else if(rank%2==1 && rank!=D->N-1)
       MPI_Send(D->plusZC,D->numPlusZC, MPI_DOUBLE, myrank+D->M, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
}


void MPI_TransferF_DSX_Zminus(Domain *D,
       double ***f1,double ***f2,double ***f3,
       double ***f4,double ***f5,double ***f6,int nx,int ny,int share)
{
    int i,j,k,numberData,start,end,nySub,nzSub,ibegin,ibottom;
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
    //first 3 is 3 layer, second 6 is 6 field variables
//    share=3;
//    numberData=share*6*(nx+5)*(nySub+5); 
    rank=(int)(myrank/D->M);   

    //Transferring even ~ odd cores 
    start=0; 
    for(k=0; k<share; k++)	// 3 layers in z direction
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f1[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f2[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f3[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f4[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f5[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f6[i][j][kstart+k];
        start+=nx;
      }

    if(rank%2==0 && rank!=D->N-1)
    {
      MPI_Recv(D->minusZ,D->numMinusZ, MPI_DOUBLE, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);  
      start=0; 
      for(k=0; k<share; k++)	// 3 layers in z direction
        for(j=0; j<ny; j++)
        {
          for(i=ibegin; i<ibottom; i++)
            f1[i][j][kend+k]=D->minusZ[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f2[i][j][kend+k]=D->minusZ[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f3[i][j][kend+k]=D->minusZ[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f4[i][j][kend+k]=D->minusZ[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f5[i][j][kend+k]=D->minusZ[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f6[i][j][kend+k]=D->minusZ[start+i];
          start+=nx;
        }  
    }
    else if(rank%2==1)
       MPI_Send(D->minusZ,D->numMinusZ, MPI_DOUBLE, myrank-D->M, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(k=0; k<share; k++)	// 3 layers in z direction
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f1[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f2[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f3[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f4[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f5[i][j][kstart+k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f6[i][j][kstart+k];
        start+=nx;
      }
        
    if(rank%2==1 && rank!=D->N-1)
    {
      MPI_Recv(D->minusZ,D->numMinusZ, MPI_DOUBLE, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);  
      start=0; 
      for(k=0; k<share; k++)	// 3 layers in z direction
        for(j=0; j<ny; j++)
        {
          for(i=ibegin; i<ibottom; i++)
            f1[i][j][kend+k]=D->minusZ[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f2[i][j][kend+k]=D->minusZ[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f3[i][j][kend+k]=D->minusZ[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f4[i][j][kend+k]=D->minusZ[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f5[i][j][kend+k]=D->minusZ[start+i];
          start+=nx;
          for(i=ibegin; i<ibottom; i++)
            f6[i][j][kend+k]=D->minusZ[start+i];
          start+=nx;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(D->minusZ,D->numMinusZ, MPI_DOUBLE, myrank-D->M, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

}


void MPI_TransferF_DSX_Zplus(Domain *D
          ,double ***f1,double ***f2,double ***f3
          ,double ***f4,double ***f5,double ***f6,int nx,int ny,int share)
{
    int i,j,k,numberData,start,end,nySub,nzSub,ibegin,ibottom;
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
    //first 2 is 2 layer, 2nd 6 is 6 field variables
//    share=3;
//    numberData=(share-1)*6*(nx+5)*(nySub+5); 
    rank=(int)(myrank/D->M);   

    //Transferring even ~ odd cores 
    start=0; 
    for(k=1; k<share; k++)
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f1[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f2[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f3[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f4[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f5[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f6[i][j][kend-k];
        start+=nx;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(D->plusZ,D->numPlusZ, MPI_DOUBLE, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);  
       start=0;
       for(k=1; k<share; k++)
         for(j=0; j<ny; j++)
         {
           for(i=ibegin; i<ibottom; i++)
             f1[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f2[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f3[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f4[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f5[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f6[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx;
         }
    }
    else if(rank%2==0 && rank!=D->N-1)
       MPI_Send(D->plusZ,D->numPlusZ, MPI_DOUBLE, myrank+D->M, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(k=1; k<share; k++)
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f1[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f2[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f3[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f4[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f5[i][j][kend-k];
        start+=nx;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f6[i][j][kend-k];
        start+=nx;
      }
        
    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(D->plusZ,D->numPlusZ, MPI_DOUBLE, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);  
       start=0;
       for(k=1; k<share; k++)
         for(j=0; j<ny; j++)
         {
           for(i=ibegin; i<ibottom; i++)
             f1[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f2[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f3[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f4[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f5[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx;
           for(i=ibegin; i<ibottom; i++)
             f6[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx;
         }
    }
    else if(rank%2==1 && rank!=D->N-1)
       MPI_Send(D->plusZ,D->numPlusZ, MPI_DOUBLE, myrank+D->M, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
}

