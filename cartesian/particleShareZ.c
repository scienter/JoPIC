#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>

void MPI_TransferP_Zplus(Domain *D,Particle ***particle,int nSpecies);
void MPI_TransferP_Zminus(Domain *D,Particle ***particle,int nSpecies);
void MPI_TransferP_Period_Z(Domain *D,Particle ***particle,int nSpecies);

void particleShareZ(Domain D)
{
  MPI_TransferP_Zplus(&D,D.particle,D.nSpecies);
  MPI_TransferP_Zminus(&D,D.particle,D.nSpecies);
  MPI_TransferP_Zplus(&D,D.track,1);
  MPI_TransferP_Zminus(&D,D.track,1);

  if(D.Period==ON && D.dimension>2)  {
    MPI_TransferP_Period_Z(&D,D.particle,D.nSpecies);
    MPI_TransferP_Period_Z(&D,D.track,1);
  }  else ;
}

void MPI_TransferP_Period_Z(Domain *D,Particle ***particle,int nSpecies)
{
  int i,j,k,n,s,numP,cnt,totalData,sendData=20,nxSub,nySub,nzSub;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks, rank;    
  double *data;
  ptclList *p,*tmp,*New;
  MPI_Status status;         

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  nxSub=D->nxSub;       nySub=0;            nzSub=0;    
  istart=D->istart;     iend=D->iend;
  jstart=0;             jend=1;
  kstart=0;             kend=1;  
	if(D->dimension>1) {	jstart=D->jstart;   jend=D->jend;   nySub=D->nySub; } else ;
	if(D->dimension>2) {  kstart=D->kstart;   kend=D->kend;   nzSub=D->nzSub; } else ;

   rank=(myrank%(D->M*D->N))/D->M;

  //Front -> Back
  numP=0;
  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      for(s=0; s<nSpecies; s++)
      {
        p=particle[i][j][kend].head[s]->pt;
        while(p)   {
          p=p->next;
          numP++;
        } 
      }

  if(D->N>1) {
    if(rank==D->N-1)
      MPI_Send(&numP,1, MPI_INT,myrank-(D->N-1)*D->M,myrank, MPI_COMM_WORLD);    
    else if(rank==0) 
      MPI_Recv(&numP,1, MPI_INT,myrank+(D->N-1)*D->M,myrank+(D->N-1)*D->M, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    
  } else ;

  totalData=numP*sendData;
  data=(double *)malloc(totalData*sizeof(double ));             

  if(D->N==1)
  {    
    n=0;
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<nSpecies; s++)
        {
          p=particle[i][j][kend].head[s]->pt;
          while(p)   
          {
            data[n*sendData+0]=p->x+i;     
            data[n*sendData+1]=p->oldX;  
            data[n*sendData+2]=p->y+j;  
            data[n*sendData+3]=p->oldY;  
            data[n*sendData+4]=p->z;  
            data[n*sendData+5]=p->oldZ-nzSub;  
            data[n*sendData+6]=p->p1;  
            data[n*sendData+7]=p->p2;  
            data[n*sendData+8]=p->p3;  
            data[n*sendData+9]=p->index;  
            data[n*sendData+10]=(double)s;  
            data[n*sendData+11]=(double)(p->core);  
            data[n*sendData+12]=p->weight;  
            data[n*sendData+13]=p->charge;  
            data[n*sendData+14]=p->p1Old1;  
            data[n*sendData+15]=p->p2Old1;  
            data[n*sendData+16]=p->p3Old1;  
            data[n*sendData+17]=p->oldX2;  
            data[n*sendData+18]=p->oldY2;  
            data[n*sendData+19]=p->oldZ2-nzSub;  
            p=p->next;
            n++;
          }
        }
    for(n=0; n<numP; n++)
    {
      i=(int)(data[n*sendData+0]);
      j=(int)(data[n*sendData+2]);
      s=(int)(data[n*sendData+10]);
      New = (ptclList *)malloc(sizeof(ptclList)); 
      New->next = particle[i][j][kstart].head[s]->pt;
      particle[i][j][kstart].head[s]->pt = New;             
      New->x=data[n*sendData+0]-i;     
      New->oldX=data[n*sendData+1];     
      New->y=data[n*sendData+2]-j;     
      New->oldY=data[n*sendData+3];     
      New->z=data[n*sendData+4];     
      New->oldZ=data[n*sendData+5];     
      New->p1=data[n*sendData+6];
      New->p2=data[n*sendData+7];
      New->p3=data[n*sendData+8];
      New->index=data[n*sendData+9];
      New->core=(int)(data[n*sendData+11]);
      New->weight=data[n*sendData+12];
      New->charge=data[n*sendData+13];
      New->p1Old1=data[n*sendData+14];
      New->p2Old1=data[n*sendData+15];
      New->p3Old1=data[n*sendData+16];
      New->oldX2=data[n*sendData+17];
      New->oldY2=data[n*sendData+18];
      New->oldZ2=data[n*sendData+19];
    }
  } else {
    if(rank==D->N-1)  
    {
      n=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(s=0; s<nSpecies; s++)
          {
            p=particle[i][j][kend].head[s]->pt;
            while(p)   
            {
              data[n*sendData+0]=p->x+i;     
              data[n*sendData+1]=p->oldX;  
              data[n*sendData+2]=p->y+j;  
              data[n*sendData+3]=p->oldY;  
              data[n*sendData+4]=p->z;  
              data[n*sendData+5]=p->oldZ-nzSub;  
              data[n*sendData+6]=p->p1;  
              data[n*sendData+7]=p->p2;  
              data[n*sendData+8]=p->p3;  
              data[n*sendData+9]=p->index;  
              data[n*sendData+10]=(double)s;  
              data[n*sendData+11]=(double)(p->core);  
              data[n*sendData+12]=p->weight;  
              data[n*sendData+13]=p->charge;  
              data[n*sendData+14]=p->p1Old1;  
              data[n*sendData+15]=p->p2Old1;  
              data[n*sendData+16]=p->p3Old1;  
              data[n*sendData+17]=p->oldX2;  
              data[n*sendData+18]=p->oldY2;  
              data[n*sendData+19]=p->oldZ2-nzSub;  
              p=p->next;
              n++;
            }
          }  

      MPI_Send(data,totalData, MPI_DOUBLE,myrank-(D->N-1)*D->M, myrank, MPI_COMM_WORLD); 
    }
    else if(rank==0)
    {    
      MPI_Recv(data,totalData, MPI_DOUBLE,myrank+(D->N-1)*D->M,myrank+(D->N-1)*D->M, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
        i=(int)(data[n*sendData+0]);
        j=(int)(data[n*sendData+2]);
        s=(int)(data[n*sendData+10]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][j][kstart].head[s]->pt;
        particle[i][j][kstart].head[s]->pt = New;             
        New->x=data[n*sendData+0]-i;     
        New->oldX=data[n*sendData+1];     
        New->y=data[n*sendData+2]-j;     
        New->oldY=data[n*sendData+3];     
        New->z=data[n*sendData+4];     
        New->oldZ=data[n*sendData+5];     
        New->p1=data[n*sendData+6];
        New->p2=data[n*sendData+7];
        New->p3=data[n*sendData+8];
        New->index=data[n*sendData+9];
        New->core=(int)(data[n*sendData+11]);
        New->weight=data[n*sendData+12];
        New->charge=data[n*sendData+13];
        New->p1Old1=data[n*sendData+14];
        New->p2Old1=data[n*sendData+15];
        New->p3Old1=data[n*sendData+16];
        New->oldX2=data[n*sendData+17];
        New->oldY2=data[n*sendData+18];
        New->oldZ2=data[n*sendData+19];
      }
    }    
    MPI_Barrier(MPI_COMM_WORLD);
  }
  free(data);

  //Back -> Front
  numP=0;
  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      for(s=0; s<nSpecies; s++)
      { 
        p=particle[i][j][kstart-1].head[s]->pt;
        while(p)   {
          p=p->next;
          numP++;
        } 
      }
  if(D->N>1) {
    if(rank==0)
      MPI_Send(&numP,1, MPI_INT,myrank+(D->N-1)*D->M,myrank, MPI_COMM_WORLD);    
    else if(rank==0) 
      MPI_Recv(&numP,1, MPI_INT,myrank-(D->N-1)*D->M,myrank-(D->N-1)*D->M, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    
  } else ;  
    
  totalData=numP*sendData;
  data=(double *)malloc(totalData*sizeof(double ));             

	if(D->N==1) 
  {
    n=0;
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<nSpecies; s++)
        {
          p=particle[i][j][kstart-1].head[s]->pt;
          while(p)   
          {
            data[n*sendData+0]=p->x+i;     
            data[n*sendData+1]=p->oldX;  
            data[n*sendData+2]=p->y+j;  
            data[n*sendData+3]=p->oldY;  
            data[n*sendData+4]=p->z;  
            data[n*sendData+5]=p->oldZ;  
            data[n*sendData+6]=p->p1;  
            data[n*sendData+7]=p->p2;  
            data[n*sendData+8]=p->p3;  
            data[n*sendData+9]=p->index;  
            data[n*sendData+10]=(double)s;  
            data[n*sendData+11]=(double)(p->core);  
            data[n*sendData+12]=p->weight;  
            data[n*sendData+13]=p->charge;  
            data[n*sendData+14]=p->p1Old1;  
            data[n*sendData+15]=p->p2Old1;  
            data[n*sendData+16]=p->p3Old1;  
            data[n*sendData+17]=p->oldX2;  
            data[n*sendData+18]=p->oldY2;  
            data[n*sendData+19]=p->oldZ2;  
            p=p->next;
            n++;
          }
        }
    for(n=0; n<numP; n++)
    {
      i=(int)(data[n*sendData+0]);
      j=(int)(data[n*sendData+2]);
      s=(int)(data[n*sendData+10]);
      New = (ptclList *)malloc(sizeof(ptclList)); 
      New->next = particle[i][j][kend-1].head[s]->pt;
      particle[i][j][kend-1].head[s]->pt = New;             
      New->x=data[n*sendData+0]-i;     
      New->oldX=data[n*sendData+1];     
      New->y=data[n*sendData+2]-j;     
      New->oldY=data[n*sendData+3];     
      New->z=data[n*sendData+4];     
      New->oldZ=data[n*sendData+5]+nzSub;     
      New->p1=data[n*sendData+6];
      New->p2=data[n*sendData+7];
      New->p3=data[n*sendData+8];
      New->index=data[n*sendData+9];
      New->core=(int)(data[n*sendData+11]);
      New->weight=data[n*sendData+12];
      New->charge=data[n*sendData+13];
      New->p1Old1=data[n*sendData+14];
      New->p2Old1=data[n*sendData+15];
      New->p3Old1=data[n*sendData+16];
      New->oldX2=data[n*sendData+17];
      New->oldY2=data[n*sendData+18];
      New->oldZ2=data[n*sendData+19]+nzSub;
    }
  } else {
    if(rank==0)
    {
      n=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(s=0; s<nSpecies; s++)
          {
            p=particle[i][j][kstart-1].head[s]->pt;
            while(p)   
            {
              data[n*sendData+0]=p->x+i;     
              data[n*sendData+1]=p->oldX;  
              data[n*sendData+2]=p->y+j;  
              data[n*sendData+3]=p->oldY;  
              data[n*sendData+4]=p->z;  
              data[n*sendData+5]=p->oldZ;  
              data[n*sendData+6]=p->p1;  
              data[n*sendData+7]=p->p2;  
              data[n*sendData+8]=p->p3;  
              data[n*sendData+9]=p->index;  
              data[n*sendData+10]=(double)s;  
              data[n*sendData+11]=(double)(p->core);  
              data[n*sendData+12]=p->weight;  
              data[n*sendData+13]=p->charge;  
              data[n*sendData+14]=p->p1Old1;  
              data[n*sendData+15]=p->p2Old1;  
              data[n*sendData+16]=p->p3Old1;  
              data[n*sendData+17]=p->oldX2;  
              data[n*sendData+18]=p->oldY2;  
              data[n*sendData+19]=p->oldZ2;  
              p=p->next;
              n++;
            }
          }
      MPI_Send(data,totalData, MPI_DOUBLE,myrank+(D->N-1)*D->M,myrank, MPI_COMM_WORLD); 
    }
    else if(rank==D->N-1) 
    {    
      MPI_Recv(data,totalData, MPI_DOUBLE,myrank-(D->N-1)*D->M,myrank-(D->N-1)*D->M,MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
        i=(int)(data[n*sendData+0]);
        j=(int)(data[n*sendData+2]);
        s=(int)(data[n*sendData+10]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][j][kend-1].head[s]->pt;
        particle[i][j][kend-1].head[s]->pt = New;             
        New->x=data[n*sendData+0]-i;     
        New->oldX=data[n*sendData+1];     
        New->y=data[n*sendData+2]-j;     
        New->oldY=data[n*sendData+3];     
        New->z=data[n*sendData+4];     
        New->oldZ=data[n*sendData+5]+nzSub;     
        New->p1=data[n*sendData+6];
        New->p2=data[n*sendData+7];
        New->p3=data[n*sendData+8];
        New->index=data[n*sendData+9];
        New->core=(int)(data[n*sendData+11]);
        New->weight=data[n*sendData+12];
        New->charge=data[n*sendData+13];
        New->p1Old1=data[n*sendData+14];
        New->p2Old1=data[n*sendData+15];
        New->p3Old1=data[n*sendData+16];
        New->oldX2=data[n*sendData+17];
        New->oldY2=data[n*sendData+18];
        New->oldZ2=data[n*sendData+19]+nzSub;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  free(data);
}

void MPI_TransferP_Zplus(Domain *D,Particle ***particle,int nSpecies)
{
   int i,j,k,n,s,numP,cnt,totalData,sendData=20,nxSub,nySub,nzSub;
  int istart,iend,jstart,jend,kstart,kend;
   int myrank, nTasks, rank;    
   double *data;
   ptclList *p,*tmp,*New;
   MPI_Status status;         

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  nxSub=D->nxSub;       nySub=0;            nzSub=0;    
  istart=D->istart;     iend=D->iend;
  jstart=0;             jend=1;
  kstart=0;             kend=1;  
	if(D->dimension>1) {	jstart=D->jstart;   jend=D->jend;   nySub=D->nySub; } else ;
	if(D->dimension>2) {  kstart=D->kstart;   kend=D->kend;   nzSub=D->nzSub; } else ;

  rank=(myrank%(D->M*D->N))/D->M;


  //Even -> odd
  if(rank%2==0 && rank!=D->N-1)
  {
    numP=0;
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<nSpecies; s++)
        {
          p=particle[i][j][kend].head[s]->pt;
          while(p)   {
            p=p->next;
            numP++;
          } 
        }
    MPI_Send(&numP,1, MPI_INT,D->nextZrank,myrank, MPI_COMM_WORLD);    
  }    
  else if(rank%2==1) 
    MPI_Recv(&numP,1, MPI_INT,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);      
  MPI_Barrier(MPI_COMM_WORLD);    

  totalData=numP*sendData;
  data=(double *)malloc(totalData*sizeof(double ));             

  if(rank%2==0 && rank!=D->N-1)
  {    
    n=0;
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<nSpecies; s++)
        {
          p=particle[i][j][kend].head[s]->pt;
          while(p)   
          {
            data[n*sendData+0]=p->x+i;     
            data[n*sendData+1]=p->oldX;  
            data[n*sendData+2]=p->y+j;  
            data[n*sendData+3]=p->oldY;  
            data[n*sendData+4]=p->z;  
            data[n*sendData+5]=p->oldZ-nzSub;  
            data[n*sendData+6]=p->p1;  
            data[n*sendData+7]=p->p2;  
            data[n*sendData+8]=p->p3;  
            data[n*sendData+9]=p->index;  
            data[n*sendData+10]=(double)s;  
            data[n*sendData+11]=(double)(p->core);  
            data[n*sendData+12]=p->weight;  
            data[n*sendData+13]=p->charge;  
            data[n*sendData+14]=p->p1Old1;  
            data[n*sendData+15]=p->p2Old1;  
            data[n*sendData+16]=p->p3Old1;  
            data[n*sendData+17]=p->oldX2;  
            data[n*sendData+18]=p->oldY2;  
            data[n*sendData+19]=p->oldZ2-nzSub;  
            p=p->next;
            n++;
          }
        }
    MPI_Send(data,totalData, MPI_DOUBLE,D->nextZrank, myrank, MPI_COMM_WORLD); 
  }

  else if(rank%2==1) 
  {    
    MPI_Recv(data,totalData, MPI_DOUBLE,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);      
    for(n=0; n<numP; n++)
    {
      i=(int)(data[n*sendData+0]);
      j=(int)(data[n*sendData+2]);
      s=(int)(data[n*sendData+10]);
      New = (ptclList *)malloc(sizeof(ptclList)); 
      New->next = particle[i][j][kstart].head[s]->pt;
      particle[i][j][kstart].head[s]->pt = New;             
      New->x=data[n*sendData+0]-i;     
      New->oldX=data[n*sendData+1];     
      New->y=data[n*sendData+2]-j;     
      New->oldY=data[n*sendData+3];     
      New->z=data[n*sendData+4];     
      New->oldZ=data[n*sendData+5];     
      New->p1=data[n*sendData+6];
      New->p2=data[n*sendData+7];
      New->p3=data[n*sendData+8];
      New->index=data[n*sendData+9];
      New->core=(int)(data[n*sendData+11]);
      New->weight=data[n*sendData+12];
      New->charge=data[n*sendData+13];
      New->p1Old1=data[n*sendData+14];
      New->p2Old1=data[n*sendData+15];
      New->p3Old1=data[n*sendData+16];
      New->oldX2=data[n*sendData+17];
      New->oldY2=data[n*sendData+18];
      New->oldZ2=data[n*sendData+19];
    }
  }
  free(data);
  MPI_Barrier(MPI_COMM_WORLD);

  //Odd -> evem
  if(rank%2==1 && rank!=D->N-1)
  {
    numP=0;
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<nSpecies; s++)
        {
          p=particle[i][j][kend].head[s]->pt;
          while(p)   {
            p=p->next;
            numP++;
          }
        } 
    MPI_Send(&numP,1, MPI_INT,D->nextZrank, myrank, MPI_COMM_WORLD);    
  }
  else if(rank%2==0 && rank!=0) 
    MPI_Recv(&numP,1, MPI_INT,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);      
  MPI_Barrier(MPI_COMM_WORLD);    

  totalData=numP*sendData;
  data=(double *)malloc(totalData*sizeof(double ));             

  if(rank%2==1 && rank!=D->N-1)
  {
    n=0;
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<nSpecies; s++)
        {
          p=particle[i][j][kend].head[s]->pt;
          while(p)   
          {
            data[n*sendData+0]=p->x+i;     
            data[n*sendData+1]=p->oldX;  
            data[n*sendData+2]=p->y+j;  
            data[n*sendData+3]=p->oldY;  
            data[n*sendData+4]=p->z;     
            data[n*sendData+5]=p->oldZ-nzSub;  
            data[n*sendData+6]=p->p1;  
            data[n*sendData+7]=p->p2;  
            data[n*sendData+8]=p->p3;  
            data[n*sendData+9]=p->index;   
            data[n*sendData+10]=(double)s;   
            data[n*sendData+11]=(double)(p->core);   
            data[n*sendData+12]=p->weight;  
            data[n*sendData+13]=p->charge;  
            data[n*sendData+14]=p->p1Old1;  
            data[n*sendData+15]=p->p2Old1;  
            data[n*sendData+16]=p->p3Old1;  
            data[n*sendData+17]=p->oldX2;  
            data[n*sendData+18]=p->oldY2;  
            data[n*sendData+19]=p->oldZ2-nzSub;  
            p=p->next;
            n++;
          }
      }
    MPI_Send(data,totalData, MPI_DOUBLE,D->nextZrank, myrank, MPI_COMM_WORLD); 
  }
 
  else if(rank%2==0 && rank!=0) 
  {    
    MPI_Recv(data,totalData, MPI_DOUBLE,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);      
    for(n=0; n<numP; n++)
    {
      i=(int)(data[n*sendData+0]);
      j=(int)(data[n*sendData+2]);
      s=(int)(data[n*sendData+10]);
      New = (ptclList *)malloc(sizeof(ptclList)); 
      New->next = particle[i][j][kstart].head[s]->pt;
      particle[i][j][kstart].head[s]->pt = New;             
      New->x=data[n*sendData+0]-i;     
      New->oldX=data[n*sendData+1];     
      New->y=data[n*sendData+2]-j;     
      New->oldY=data[n*sendData+3];     
      New->z=data[n*sendData+4];     
      New->oldZ=data[n*sendData+5];     
      New->p1=data[n*sendData+6];
      New->p2=data[n*sendData+7];
      New->p3=data[n*sendData+8];
      New->index=data[n*sendData+9];
      New->core=(int)(data[n*sendData+11]);
      New->weight=data[n*sendData+12];
      New->charge=data[n*sendData+13];
      New->p1Old1=data[n*sendData+14];
      New->p2Old1=data[n*sendData+15];
      New->p3Old1=data[n*sendData+16];
      New->oldX2=data[n*sendData+17];
      New->oldY2=data[n*sendData+18];
      New->oldZ2=data[n*sendData+19];
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  free(data);

}


void MPI_TransferP_Zminus(Domain *D,Particle ***particle,int nSpecies)
{
   int i,j,k,n,s,numP,cnt,totalData,sendData=20,nxSub,nySub,nzSub;
  int istart,iend,jstart,jend,kstart,kend;
   int myrank, nTasks, rank;
   double *data;
   ptclList *p,*tmp,*New;
   MPI_Status status;         

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  nxSub=D->nxSub;       nySub=0;            nzSub=0;    
  istart=D->istart;     iend=D->iend;
  jstart=0;             jend=1;
  kstart=0;             kend=1;  
	if(D->dimension>1) {	jstart=D->jstart;   jend=D->jend;   nySub=D->nySub; } else ;
	if(D->dimension>2) {  kstart=D->kstart;   kend=D->kend;   nzSub=D->nzSub; } else ;

   rank=(myrank%(D->M*D->N))/D->M;

    //Even -> odd
    if(rank%2==0 && rank!=0)
    {
      numP=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(s=0; s<nSpecies; s++)
          { 
            p=particle[i][j][kstart-1].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            } 
          }
      MPI_Send(&numP,1, MPI_INT,D->prevZrank, myrank, MPI_COMM_WORLD);    
    }    
    else if(rank%2==1 && rank!=D->N-1) 
      MPI_Recv(&numP,1, MPI_INT,D->nextZrank,D->nextZrank,MPI_COMM_WORLD,&status);        
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    data=(double *)malloc(totalData*sizeof(double ));             
 
    if(rank%2==0 && rank!=0)
    {
      n=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(s=0; s<nSpecies; s++)
          {
            p=particle[i][j][kstart-1].head[s]->pt;
            while(p)   
            {
              data[n*sendData+0]=p->x+i;     
              data[n*sendData+1]=p->oldX;  
              data[n*sendData+2]=p->y+j;  
              data[n*sendData+3]=p->oldY;  
              data[n*sendData+4]=p->z;  
              data[n*sendData+5]=p->oldZ;  
              data[n*sendData+6]=p->p1;  
              data[n*sendData+7]=p->p2;  
              data[n*sendData+8]=p->p3;  
              data[n*sendData+9]=p->index;  
              data[n*sendData+10]=(double)s;  
              data[n*sendData+11]=(double)(p->core);  
              data[n*sendData+12]=p->weight;  
              data[n*sendData+13]=p->charge;  
              data[n*sendData+14]=p->p1Old1;  
              data[n*sendData+15]=p->p2Old1;  
              data[n*sendData+16]=p->p3Old1;  
              data[n*sendData+17]=p->oldX2;  
              data[n*sendData+18]=p->oldY2;  
              data[n*sendData+19]=p->oldZ2;  
              p=p->next;
              n++;
            }
          }
      MPI_Send(data,totalData, MPI_DOUBLE,D->prevZrank,myrank, MPI_COMM_WORLD); 
    }
    else if(rank%2==1 && rank!=D->N-1) 
    {    
      MPI_Recv(data,totalData, MPI_DOUBLE,D->nextZrank,D->nextZrank,MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
        i=(int)(data[n*sendData+0]);
        j=(int)(data[n*sendData+2]);
        s=(int)(data[n*sendData+10]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][j][kend-1].head[s]->pt;
        particle[i][j][kend-1].head[s]->pt = New;             
        New->x=data[n*sendData+0]-i;     
        New->oldX=data[n*sendData+1];     
        New->y=data[n*sendData+2]-j;     
        New->oldY=data[n*sendData+3];     
        New->z=data[n*sendData+4];     
        New->oldZ=data[n*sendData+5]+nzSub;     
        New->p1=data[n*sendData+6];
        New->p2=data[n*sendData+7];
        New->p3=data[n*sendData+8];
        New->index=data[n*sendData+9];
        New->core=(int)(data[n*sendData+11]);
        New->weight=data[n*sendData+12];
        New->charge=data[n*sendData+13];
        New->p1Old1=data[n*sendData+14];
        New->p2Old1=data[n*sendData+15];
        New->p3Old1=data[n*sendData+16];
        New->oldX2=data[n*sendData+17];
        New->oldY2=data[n*sendData+18];
        New->oldZ2=data[n*sendData+19]+nzSub;
      }
    }
    free(data);
    MPI_Barrier(MPI_COMM_WORLD);

    //Odd -> evem
    if(rank%2==1)
    {
      numP=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(s=0; s<nSpecies; s++)
          {
            p=particle[i][j][kstart-1].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            }
          } 
      MPI_Send(&numP,1, MPI_INT,D->prevZrank, myrank, MPI_COMM_WORLD);    
    }
    else if(rank%2==0 && rank!=D->N-1) 
      MPI_Recv(&numP,1, MPI_INT,D->nextZrank,D->nextZrank,MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    data=(double *)malloc(totalData*sizeof(double ));             

    if(rank%2==1)
    {
      n=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(s=0; s<nSpecies; s++)
          {
            p=particle[i][j][kstart-1].head[s]->pt;
            while(p)   
            {
              data[n*sendData+0]=p->x+i;     
              data[n*sendData+1]=p->oldX;  
              data[n*sendData+2]=p->y+j;  
              data[n*sendData+3]=p->oldY;  
              data[n*sendData+4]=p->z;  
              data[n*sendData+5]=p->oldZ;  
              data[n*sendData+6]=p->p1;  
              data[n*sendData+7]=p->p2;  
              data[n*sendData+8]=p->p3;  
              data[n*sendData+9]=p->index;  
              data[n*sendData+10]=(double)s;  
              data[n*sendData+11]=(double)(p->core);  
              data[n*sendData+12]=p->weight;  
              data[n*sendData+13]=p->charge;  
              data[n*sendData+14]=p->p1Old1;  
              data[n*sendData+15]=p->p2Old1;  
              data[n*sendData+16]=p->p3Old1;  
              data[n*sendData+17]=p->oldX2;  
              data[n*sendData+18]=p->oldY2;  
              data[n*sendData+19]=p->oldZ2;  
              p=p->next;
              n++;
            }
          }
      MPI_Send(data,totalData, MPI_DOUBLE,D->prevZrank, myrank, MPI_COMM_WORLD); 
    }

    else if(rank%2==0 && rank!=D->N-1) 
    {    
      MPI_Recv(data,totalData, MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
        i=(int)(data[n*sendData+0]);
        j=(int)(data[n*sendData+2]);
        s=(int)(data[n*sendData+10]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][j][kend-1].head[s]->pt;
        particle[i][j][kend-1].head[s]->pt = New;             
        New->x=data[n*sendData+0]-i;     
        New->oldX=data[n*sendData+1];     
        New->y=data[n*sendData+2]-j;     
        New->oldY=data[n*sendData+3];     
        New->z=data[n*sendData+4];     
        New->oldZ=data[n*sendData+5]+nzSub;     
        New->p1=data[n*sendData+6];
        New->p2=data[n*sendData+7];
        New->p3=data[n*sendData+8];
        New->index=data[n*sendData+9];
        New->core=(int)(data[n*sendData+11]);
        New->weight=data[n*sendData+12];
        New->charge=data[n*sendData+13];
        New->p1Old1=data[n*sendData+14];
        New->p2Old1=data[n*sendData+15];
        New->p3Old1=data[n*sendData+16];
        New->oldX2=data[n*sendData+17];
        New->oldY2=data[n*sendData+18];
        New->oldZ2=data[n*sendData+19]+nzSub;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);

}

