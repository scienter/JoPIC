#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>

void particleShareZ(Domain *D)
{
  int istart,iend,jstart,jend,kstart,kend;
  void MPI_TransferP_Zplus();
  void MPI_TransferP_Zminus();

  istart=D->istart-1;
  iend=D->iend+1;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;

  MPI_TransferP_Zplus(D,istart,iend,jstart,jend,kstart,kend);
  MPI_TransferP_Zminus(D,istart,iend,jstart,jend,kstart,kend);
}

void MPI_TransferP_Zplus(Domain *D
    ,int istart,int iend,int jstart,int jend,int kstart,int kend)
{
   int i,j,k,n,s,numP,cnt,totalData,sendData=12,nxSub,nySub,nzSub;
   int myrank, nTasks, rank;
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

   rank=(int)(myrank/D->M);

    //Even -> odd
    if(rank%2==0 && rank!=D->N-1)
    {
      numP=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][kend].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            } 
          }
      MPI_Send(&numP,1, MPI_INT, myrank+D->M, myrank, MPI_COMM_WORLD);    
    }    
    else if(rank%2==1) 
      MPI_Recv(&numP,1, MPI_INT, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    upP=(double *)malloc(totalData*sizeof(double ));             

    if(rank%2==0 && rank!=D->N-1)
    {    
      n=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][kend].head[s]->pt;
            while(p)   
            {
              upP[n*sendData+0]=p->x+i;     
              upP[n*sendData+1]=p->oldX;  
              upP[n*sendData+2]=p->y+j;  
              upP[n*sendData+3]=p->oldY;  
              upP[n*sendData+4]=p->z;  
              upP[n*sendData+5]=p->oldZ-nzSub;  
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
      MPI_Send(upP,totalData, MPI_DOUBLE, myrank+D->M, myrank, MPI_COMM_WORLD); 
    }

    else if(rank%2==1) 
    {    
      MPI_Recv(upP,totalData, MPI_DOUBLE, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
        i=(int)(upP[n*sendData+0]);
        j=(int)(upP[n*sendData+2]);
        s=(int)(upP[n*sendData+10]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][j][kstart].head[s]->pt;
        particle[i][j][kstart].head[s]->pt = New;             
        New->x=upP[n*sendData+0]-i;     
        New->oldX=upP[n*sendData+1];     
        New->y=upP[n*sendData+2]-j;     
        New->oldY=upP[n*sendData+3];     
        New->z=upP[n*sendData+4];     
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

    //Odd -> evem
    if(rank%2==1 && rank!=D->N-1)
    {
      numP=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][kend].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            }
          } 
      MPI_Send(&numP,1, MPI_INT, myrank+D->M, myrank, MPI_COMM_WORLD);    
    }
    else if(rank%2==0 && rank!=0) 
      MPI_Recv(&numP,1, MPI_INT, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    upP=(double *)malloc(totalData*sizeof(double ));             

    if(rank%2==1 && rank!=D->N-1)
    {
      n=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][kend].head[s]->pt;
            while(p)   
            {
              upP[n*sendData+0]=p->x+i;     
              upP[n*sendData+1]=p->oldX;  
              upP[n*sendData+2]=p->y+j;  
              upP[n*sendData+3]=p->oldY;  
              upP[n*sendData+4]=p->z;     
              upP[n*sendData+5]=p->oldZ-nzSub;  
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
      MPI_Send(upP,totalData, MPI_DOUBLE, myrank+D->M, myrank, MPI_COMM_WORLD); 
    }
 
    else if(rank%2==0 && rank!=0) 
    {    
      MPI_Recv(upP,totalData, MPI_DOUBLE, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
        i=(int)(upP[n*sendData+0]);
        j=(int)(upP[n*sendData+2]);
        s=(int)(upP[n*sendData+10]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][j][kstart].head[s]->pt;
        particle[i][j][kstart].head[s]->pt = New;             
        New->x=upP[n*sendData+0]-i;     
        New->oldX=upP[n*sendData+1];     
        New->y=upP[n*sendData+2]-j;     
        New->oldY=upP[n*sendData+3];     
        New->z=upP[n*sendData+4];     
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


void MPI_TransferP_Zminus(Domain *D
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

   rank=(int)(myrank/D->M);

    //Even -> odd
    if(rank%2==0 && rank!=0)
    {
      numP=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          { 
            p=particle[i][j][kstart-1].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            } 
          }
      MPI_Send(&numP,1, MPI_INT, myrank-D->M, myrank, MPI_COMM_WORLD);    
    }    
    else if(rank%2==1 && rank!=D->N-1) 
      MPI_Recv(&numP,1, MPI_INT, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);        
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    btP=(double *)malloc(totalData*sizeof(double ));             
 
    if(rank%2==0 && rank!=0)
    {
      n=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][kstart-1].head[s]->pt;
            while(p)   
            {
              btP[n*sendData+0]=p->x+i;     
              btP[n*sendData+1]=p->oldX;  
              btP[n*sendData+2]=p->y+j;  
              btP[n*sendData+3]=p->oldY;  
              btP[n*sendData+4]=p->z;  
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
      MPI_Send(btP,totalData, MPI_DOUBLE, myrank-D->M, myrank, MPI_COMM_WORLD); 
    }
    else if(rank%2==1 && rank!=D->N-1) 
    {    
      MPI_Recv(btP,totalData, MPI_DOUBLE, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);      
      n=0;
      while(numP>0)
      {
        i=(int)(btP[n*sendData+0]);
        j=(int)(btP[n*sendData+2]);
        s=(int)(btP[n*sendData+10]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][j][kend-1].head[s]->pt;
        particle[i][j][kend-1].head[s]->pt = New;             
        New->x=btP[n*sendData+0]-i;     
        New->oldX=btP[n*sendData+1];     
        New->y=btP[n*sendData+2]-j;     
        New->oldY=btP[n*sendData+3];     
        New->z=btP[n*sendData+4];     
        New->oldZ=btP[n*sendData+5]+nzSub;     
        New->p1=btP[n*sendData+6];
        New->p2=btP[n*sendData+7];
        New->p3=btP[n*sendData+8];
        New->index=(int)(btP[n*sendData+9]);
        New->core=(int)(btP[n*sendData+11]);
        n++;
        numP--;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(btP);

    //Odd -> evem
    if(rank%2==1)
    {
      numP=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][kstart-1].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            }
          } 
      MPI_Send(&numP,1, MPI_INT, myrank-D->M, myrank, MPI_COMM_WORLD);    
    }
    else if(rank%2==0 && rank!=D->N-1) 
      MPI_Recv(&numP,1, MPI_INT, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);      
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
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][kstart-1].head[s]->pt;
            while(p)   
            {
              btP[n*sendData+0]=p->x+i;     
              btP[n*sendData+1]=p->oldX;  
              btP[n*sendData+2]=p->y+j;  
              btP[n*sendData+3]=p->oldY;  
              btP[n*sendData+4]=p->z;  
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
      MPI_Send(btP,totalData, MPI_DOUBLE, myrank-D->M, myrank, MPI_COMM_WORLD); 
    }

    else if(rank%2==0 && rank!=D->N-1) 
    {    
      MPI_Recv(btP,totalData, MPI_DOUBLE, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);      
      n=0;
      while(numP>0)
      {
        i=(int)(btP[n*sendData+0]);
        j=(int)(btP[n*sendData+2]);
        s=(int)(btP[n*sendData+10]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][j][kend-1].head[s]->pt;
        particle[i][j][kend-1].head[s]->pt = New;             
        New->x=btP[n*sendData+0]-i;     
        New->oldX=btP[n*sendData+1];     
        New->y=btP[n*sendData+2]-j;     
        New->oldY=btP[n*sendData+3];     
        New->z=btP[n*sendData+4];     
        New->oldZ=btP[n*sendData+5]+nzSub;     
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

