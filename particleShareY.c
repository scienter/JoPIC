#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>

void particleShareY(Domain *D)
{
  int istart,iend,jstart,jend,kstart,kend;
  void MPI_TransferP_Yplus();
  void MPI_TransferP_Yminus();

  switch(D->dimension)  {
  case 1 :
    break;
  case 2 :
    istart=D->istart-1;
    iend=D->iend+1;
    jstart=D->jstart;
    jend=D->jend;
    kstart=0;
    kend=1;
    MPI_TransferP_Yplus(D,istart,iend,jstart,jend,kstart,kend);
    MPI_TransferP_Yminus(D,istart,iend,jstart,jend,kstart,kend);
    break;
  case 3 :
    MPI_TransferP_Yplus(D,D->istart-1,D->iend+1,D->jstart,D->jend,D->kstart-1,D->kend+1);
    MPI_TransferP_Yminus(D,D->istart-1,D->iend+1,D->jstart,D->jend,D->kstart-1,D->kend+1);
    break;
  default :
    printf("In particleShareY, what dimenstion(%d)?\n",D->dimension);
  }
}

void MPI_TransferP_Yplus(Domain *D
          ,int istart,int iend,int jstart,int jend,int kstart,int kend)
{
   int i,j,k,n,s,numP,cnt,totalData,sendData=12,nxSub,nySub,nzSub;   int myrank, nTasks, rank;
   Particle ***particle;
   particle=D->particle;     
   double *upP;
   ptclList *p,*tmp,*New;
   MPI_Status status;         

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

   nxSub=D->nxSub;
   nySub=D->nySub;
   nzSub=D->nzSub;

   rank=myrank%D->M;

    //Even -> odd
    if(rank%2==0 && rank!=D->M-1)
    {
      numP=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][jend][k].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            } 
          }
      MPI_Send(&numP,1, MPI_INT, myrank+1, myrank, MPI_COMM_WORLD);    
    }    
    else if(rank%2==1) 
      MPI_Recv(&numP,1, MPI_INT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    upP=(double *)malloc(totalData*sizeof(double ));             

    if(rank%2==0 && rank!=D->M-1)
    {    
      n=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][jend][k].head[s]->pt;
            while(p)   
            {
              upP[n*sendData+0]=p->x+i;     
              upP[n*sendData+1]=p->oldX;  
              upP[n*sendData+2]=p->y;  
              upP[n*sendData+3]=p->oldY-nySub;  
              upP[n*sendData+4]=p->z+k;  
              upP[n*sendData+5]=p->oldZ;  
              upP[n*sendData+6]=p->p1;  
              upP[n*sendData+7]=p->p2;  
              upP[n*sendData+8]=p->p3;  
              upP[n*sendData+9]=(double)(p->index);  
              upP[n*sendData+10]=(double)s;  
              upP[n*sendData+11]=(double)(p->core);  
              p=p->next;
              n++;
            }
          }
      MPI_Send(upP,totalData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD); 
    }

    else if(rank%2==1) 
    {    
      MPI_Recv(upP,totalData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
        i=(int)(upP[n*sendData+0]);
        k=(int)(upP[n*sendData+4]);
        s=(int)(upP[n*sendData+10]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][jstart][k].head[s]->pt;
        particle[i][jstart][k].head[s]->pt = New;             
        New->x=upP[n*sendData+0]-i;     
        New->oldX=upP[n*sendData+1];     
        New->y=upP[n*sendData+2];     
        New->oldY=upP[n*sendData+3];     
        New->z=upP[n*sendData+4]-k;     
        New->oldZ=upP[n*sendData+5];     
        New->p1=upP[n*sendData+6];
        New->p2=upP[n*sendData+7];
        New->p3=upP[n*sendData+8];
        New->index=(int)(upP[n*sendData+9]);
        New->core=(int)(upP[n*sendData+11]);
      }
    }
    free(upP);
    MPI_Barrier(MPI_COMM_WORLD);

    //Odd -> evem
    if(rank%2==1 && rank!=D->M-1)
    {
      numP=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][jend][k].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            }
          } 
      MPI_Send(&numP,1, MPI_INT, myrank+1, myrank, MPI_COMM_WORLD);    
    }
    else if(rank%2==0 && rank!=0) 
      MPI_Recv(&numP,1, MPI_INT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    upP=(double *)malloc(totalData*sizeof(double ));             

    if(rank%2==1 && rank!=D->M-1)
    {
      n=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][jend][k].head[s]->pt;
            while(p)   
            {
              upP[n*sendData+0]=p->x+i;     
              upP[n*sendData+1]=p->oldX;  
              upP[n*sendData+2]=p->y;  
              upP[n*sendData+3]=p->oldY-nySub;  
              upP[n*sendData+4]=p->z+k;     
              upP[n*sendData+5]=p->oldZ;  
              upP[n*sendData+6]=p->p1;  
              upP[n*sendData+7]=p->p2;  
              upP[n*sendData+8]=p->p3;  
              upP[n*sendData+9]=(double)(p->index);   
              upP[n*sendData+10]=(double)s;   
              upP[n*sendData+11]=(double)(p->core);   
              p=p->next;
              n++;
            }
        }
      MPI_Send(upP,totalData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD); 
    }
 
    else if(rank%2==0 && rank!=0) 
    {    
      MPI_Recv(upP,totalData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
        i=(int)(upP[n*sendData+0]);
        k=(int)(upP[n*sendData+4]);
        s=(int)(upP[n*sendData+10]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][jstart][k].head[s]->pt;
        particle[i][jstart][k].head[s]->pt = New;             
        New->x=upP[n*sendData+0]-i;     
        New->oldX=upP[n*sendData+1];     
        New->y=upP[n*sendData+2];     
        New->oldY=upP[n*sendData+3];     
        New->z=upP[n*sendData+4]-k;     
        New->oldZ=upP[n*sendData+5];     
        New->p1=upP[n*sendData+6];
        New->p2=upP[n*sendData+7];
        New->p3=upP[n*sendData+8];
        New->index=(int)(upP[n*sendData+9]);
        New->core=(int)(upP[n*sendData+11]);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(upP);
}


void MPI_TransferP_Yminus(Domain *D
          ,int istart,int iend,int jstart,int jend,int kstart,int kend)
{
   int i,j,k,n,s,numP,cnt,totalData,sendData=12,nxSub,nySub,nzSub;
   int myrank, nTasks, rank;
   Particle ***particle;
   particle=D->particle;     
   double *btP;
   ptclList *p,*tmp,*New;
   MPI_Status status;         

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

   nxSub=D->nxSub;
   nySub=D->nySub;
   nzSub=D->nzSub;

   rank=myrank%D->M;

    //Even -> odd
    if(rank%2==0 && rank!=0)
    {
      numP=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(k=kstart; k<kend; k++)
          { 
            p=particle[i][jstart-1][k].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            } 
          }
      MPI_Send(&numP,1, MPI_INT, myrank-1, myrank, MPI_COMM_WORLD);    
    }    
    else if(rank%2==1 && rank!=D->M-1) 
      MPI_Recv(&numP,1, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);        
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    btP=(double *)malloc(totalData*sizeof(double ));             
//     for(i=0; i<totalData; i++)
//       upP[i]=0.0;
 
    if(rank%2==0 && rank!=0)
    {
      n=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][jstart-1][k].head[s]->pt;
            while(p)   
            {
              btP[n*sendData+0]=p->x+i;     
              btP[n*sendData+1]=p->oldX;  
              btP[n*sendData+2]=p->y;  
              btP[n*sendData+3]=p->oldY;  
              btP[n*sendData+4]=p->z+k;  
              btP[n*sendData+5]=p->oldZ;  
              btP[n*sendData+6]=p->p1;  
              btP[n*sendData+7]=p->p2;  
              btP[n*sendData+8]=p->p3;  
              btP[n*sendData+9]=(double)(p->index);  
              btP[n*sendData+10]=(double)s;  
              btP[n*sendData+11]=(double)(p->core);  
              p=p->next;
              n++;
            }
          }
      MPI_Send(btP,totalData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD); 
    }
    else if(rank%2==1 && rank!=D->M-1) 
    {    
      MPI_Recv(btP,totalData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
      n=0;
      while(numP>0)
      {
        i=(int)(btP[n*sendData+0]);
        k=(int)(btP[n*sendData+4]);
        s=(int)(btP[n*sendData+10]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][jend-1][k].head[s]->pt;
        particle[i][jend-1][k].head[s]->pt = New;             
        New->x=btP[n*sendData+0]-i;     
        New->oldX=btP[n*sendData+1];     
        New->y=btP[n*sendData+2];     
        New->oldY=btP[n*sendData+3]+nySub;     
        New->z=btP[n*sendData+4]-k;     
        New->oldZ=btP[n*sendData+5];     
        New->p1=btP[n*sendData+6];
        New->p2=btP[n*sendData+7];
        New->p3=btP[n*sendData+8];
        New->index=(int)(btP[n*sendData+9]);
        New->core=(int)(btP[n*sendData+11]);
        n++;
        numP--;
      }
    }
    free(btP);
    MPI_Barrier(MPI_COMM_WORLD);

    //Odd -> evem
    if(rank%2==1)
    {
      numP=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][jstart-1][k].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            }
          } 
      MPI_Send(&numP,1, MPI_INT, myrank-1, myrank, MPI_COMM_WORLD);    
    }
    else if(rank%2==0 && rank!=D->M-1) 
      MPI_Recv(&numP,1, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    btP=(double *)malloc(totalData*sizeof(double ));             
//    for(i=0; i<totalData; i++)
//      btP[i]=0.0;

    if(rank%2==1)
    {
      n=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][jstart-1][k].head[s]->pt;
            while(p)   
            {
              btP[n*sendData+0]=p->x+i;     
              btP[n*sendData+1]=p->oldX;  
              btP[n*sendData+2]=p->y;  
              btP[n*sendData+3]=p->oldY;  
              btP[n*sendData+4]=p->z+k;  
              btP[n*sendData+5]=p->oldZ;  
              btP[n*sendData+6]=p->p1;  
              btP[n*sendData+7]=p->p2;  
              btP[n*sendData+8]=p->p3;  
              btP[n*sendData+9]=(double)(p->index);  
              btP[n*sendData+10]=(double)s;  
              btP[n*sendData+11]=(double)(p->core);  
              p=p->next;
              n++;
            }
          }
      MPI_Send(btP,totalData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD); 
    }

    else if(rank%2==0 && rank!=D->M-1) 
    {    
      MPI_Recv(btP,totalData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
      n=0;
      while(numP>0)
      {
        i=(int)(btP[n*sendData+0]);
        k=(int)(btP[n*sendData+4]);
        s=(int)(btP[n*sendData+10]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][jend-1][k].head[s]->pt;
        particle[i][jend-1][k].head[s]->pt = New;             
        New->x=btP[n*sendData+0]-i;     
        New->oldX=btP[n*sendData+1];     
        New->y=btP[n*sendData+2];     
        New->oldY=btP[n*sendData+3]+nySub;     
        New->z=btP[n*sendData+4]-k;     
        New->oldZ=btP[n*sendData+5];     
        New->p1=btP[n*sendData+6];
        New->p2=btP[n*sendData+7];
        New->p3=btP[n*sendData+8];
        New->index=(int)(btP[n*sendData+9]);
        New->core=(int)(btP[n*sendData+11]);
        numP--;
        n++;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(btP);

}

