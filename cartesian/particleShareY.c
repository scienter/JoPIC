#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>

void MPI_TransferP_Period_Y(Domain *D,Particle ***particle,int nSpecies);
void MPI_TransferP_Yplus(Domain *D,Particle ***particle,int nSpecies);
void MPI_TransferP_Yminus(Domain *D,Particle ***particle,int nSpecies);

void particleShareY(Domain D)
{
  MPI_TransferP_Yplus(&D,D.particle,D.nSpecies);
  MPI_TransferP_Yminus(&D,D.particle,D.nSpecies);
  MPI_TransferP_Yplus(&D,D.track,1);
  MPI_TransferP_Yminus(&D,D.track,1);
  if(D.Period==ON && D.dimension>1) {
    MPI_TransferP_Period_Y(&D,D.particle,D.nSpecies);
    MPI_TransferP_Period_Y(&D,D.track,1);
	} else ;
}

void MPI_TransferP_Period_Y(Domain *D,Particle ***particle,int nSpecies)
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
	if(D->dimension>2) {  kstart=D->kstart-1; kend=D->kend+1; nzSub=D->nzSub; } else ;

  rank=(myrank%(D->M*D->N))%D->M;

	// Up -> Bottom
  numP=0;
  for(i=istart; i<iend; i++)
   	for(k=kstart; k<kend; k++)
      for(s=0; s<nSpecies; s++)
      {
        p=particle[i][jend][k].head[s]->pt;
          while(p)   {
            p=p->next;
            numP++;
          } 
			}
	if(D->M>1) {
		if(rank==D->M-1)
			MPI_Send(&numP,1, MPI_INT,myrank-(D->M-1),myrank, MPI_COMM_WORLD);    
		else if(rank==0) 
      	MPI_Recv(&numP,1, MPI_INT,myrank+(D->M-1),myrank+(D->M-1), MPI_COMM_WORLD,&status);      
    	MPI_Barrier(MPI_COMM_WORLD);    
	} else ;

	totalData=numP*sendData;
  data=(double *)malloc(totalData*sizeof(double ));             

	if(D->M==1) {
		n=0;
   	for(i=istart; i<iend; i++)
   		for(k=kstart; k<kend; k++)
      	for(s=0; s<nSpecies; s++)
	      {	
         	p=particle[i][jend][k].head[s]->pt;
   	      while(p)   
      	  {
         	  data[n*sendData+0]=p->x+i;     
            data[n*sendData+1]=p->oldX;  
            data[n*sendData+2]=p->y;  
	          data[n*sendData+3]=p->oldY-nySub;  
   	        data[n*sendData+4]=p->z+k;  
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
         	  data[n*sendData+18]=p->oldY2-nySub;  
            data[n*sendData+19]=p->oldZ2;  
            p=p->next;
	          n++;
   	      }
      	}
		for(n=0; n<numP; n++)
    {
      i=(int)(data[n*sendData+0]);
      k=(int)(data[n*sendData+4]);
      s=(int)(data[n*sendData+10]);
      New = (ptclList *)malloc(sizeof(ptclList)); 
      New->next = particle[i][jstart][k].head[s]->pt;
      particle[i][jstart][k].head[s]->pt = New;             
      New->x=data[n*sendData+0]-i;     
      New->oldX=data[n*sendData+1];     
      New->y=data[n*sendData+2];     
      New->oldY=data[n*sendData+3];     
      New->z=data[n*sendData+4]-k;     
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
		if(rank==D->M-1) {
			n=0;
   		for(i=istart; i<iend; i++)
   			for(k=kstart; k<kend; k++)
      	  for(s=0; s<nSpecies; s++)
	        {	
          	p=particle[i][jend][k].head[s]->pt;
	   	       while(p)   
   	   	    {
      	     	data[n*sendData+0]=p->x+i;     
           	  data[n*sendData+1]=p->oldX;  
          	  data[n*sendData+2]=p->y;  
	          	data[n*sendData+3]=p->oldY-nySub;  
   	         	data[n*sendData+4]=p->z+k;  
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
             	data[n*sendData+18]=p->oldY2-nySub;  
	          	data[n*sendData+19]=p->oldZ2;  
   	         	p=p->next;
	   	        n++;
   	   	    }
      	  }

			MPI_Send(data,totalData, MPI_DOUBLE,myrank-(D->M-1),myrank, MPI_COMM_WORLD); 
		}
		else if(rank==0) 
		{    
      MPI_Recv(data,totalData, MPI_DOUBLE,myrank+(D->M-1),myrank+(D->M-1),MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
      	i=(int)(data[n*sendData+0]);
      	k=(int)(data[n*sendData+4]);
      	s=(int)(data[n*sendData+10]);
      	New = (ptclList *)malloc(sizeof(ptclList)); 
      	New->next = particle[i][jstart][k].head[s]->pt;
      	particle[i][jstart][k].head[s]->pt = New;             
      	New->x=data[n*sendData+0]-i;     
      	New->oldX=data[n*sendData+1];     
      	New->y=data[n*sendData+2];     
      	New->oldY=data[n*sendData+3];     
      	New->z=data[n*sendData+4]-k;     
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

   // Bottom -> Up
	numP=0;
   for(i=istart; i<iend; i++)
   	for(k=kstart; k<kend; k++)
      	for(s=0; s<nSpecies; s++)
         { 
         	p=particle[i][jstart-1][k].head[s]->pt;
            while(p)   {
            	p=p->next;
            	numP++;
           	} 
			}
	if(D->M>1) {
		if(rank==0)
      	MPI_Send(&numP,1, MPI_INT,myrank+(D->M-1), myrank, MPI_COMM_WORLD);    
		else if(rank==D->M-1) 
      	MPI_Recv(&numP,1, MPI_INT,myrank-(D->M-1),myrank-(D->M-1),MPI_COMM_WORLD,&status);        
		MPI_Barrier(MPI_COMM_WORLD);    
	} else ;

	totalData=numP*sendData;
	data=(double *)malloc(totalData*sizeof(double ));             
 
	if(D->M==1) 
	{
		n=0;
   	for(i=istart; i<iend; i++)
   		for(k=kstart; k<kend; k++)
      	for(s=0; s<nSpecies; s++)
        {	
        	p=particle[i][jstart-1][k].head[s]->pt;
	          while(p)   
   	       {
      	    data[n*sendData+0]=p->x+i;     
            data[n*sendData+1]=p->oldX;  
          	data[n*sendData+2]=p->y;  
            data[n*sendData+3]=p->oldY;  
	          data[n*sendData+4]=p->z+k;  
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
      k=(int)(data[n*sendData+4]);
      s=(int)(data[n*sendData+10]);
      New = (ptclList *)malloc(sizeof(ptclList)); 
      New->next = particle[i][jend-1][k].head[s]->pt;
      particle[i][jend-1][k].head[s]->pt = New;             
      New->x=data[n*sendData+0]-i;     
      New->oldX=data[n*sendData+1];     
      New->y=data[n*sendData+2];     
      New->oldY=data[n*sendData+3]+nySub;     
      New->z=data[n*sendData+4]-k;     
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
      New->oldY2=data[n*sendData+18]+nySub;
      New->oldZ2=data[n*sendData+19];
    }
	} else {
		if(rank==0)
		{
			n=0;
	   	for(i=istart; i<iend; i++)
   			for(k=kstart; k<kend; k++)
      		for(s=0; s<nSpecies; s++)
         	{	
         		p=particle[i][jstart-1][k].head[s]->pt;
	          while(p)   
	   	      {
   	   	    	data[n*sendData+0]=p->x+i;     
      	      data[n*sendData+1]=p->oldX;  
         	 	  data[n*sendData+2]=p->y;  
          	  data[n*sendData+3]=p->oldY;  
	          	data[n*sendData+4]=p->z+k;  
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
      MPI_Send(data,totalData, MPI_DOUBLE,myrank+(D->M-1),myrank, MPI_COMM_WORLD); 
		}
    else if(rank==D->M-1) 
		{    
      MPI_Recv(data,totalData, MPI_DOUBLE,myrank-(D->M-1),myrank-(D->M-1),MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
      	i=(int)(data[n*sendData+0]);
      	k=(int)(data[n*sendData+4]);
      	s=(int)(data[n*sendData+10]);
      	New = (ptclList *)malloc(sizeof(ptclList)); 
      	New->next = particle[i][jend-1][k].head[s]->pt;
      	particle[i][jend-1][k].head[s]->pt = New;             
      	New->x=data[n*sendData+0]-i;     
      	New->oldX=data[n*sendData+1];     
      	New->y=data[n*sendData+2];     
      	New->oldY=data[n*sendData+3]+nySub;     
      	New->z=data[n*sendData+4]-k;     
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
      	New->oldY2=data[n*sendData+18]+nySub;
      	New->oldZ2=data[n*sendData+19];
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
	}
   free(data);
}


void MPI_TransferP_Yplus(Domain *D,Particle ***particle,int nSpecies)
{
   int i,j,k,n,s,numP,cnt,totalData,sendData=20,nxSub,nySub,nzSub;
   int istart,iend,jstart,jend,kstart,kend;
   int myrank, nTasks, rank;   
   double *upP;
   ptclList *p,*tmp,*New;
   MPI_Status status;         

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  nxSub=D->nxSub;       nySub=0;            nzSub=0;    
  istart=D->istart;     iend=D->iend;
  jstart=0;             jend=1;
  kstart=0;             kend=1;  
	if(D->dimension>1) {	jstart=D->jstart;   jend=D->jend;   nySub=D->nySub; } else ;
	if(D->dimension>2) {  kstart=D->kstart-1; kend=D->kend+1; nzSub=D->nzSub; } else ;

   rank=(myrank%(D->M*D->N))%D->M;

    //Even -> odd
    if(rank%2==0 && rank!=D->M-1)
    {
      numP=0;
      for(i=istart; i<iend; i++)
        for(k=kstart; k<kend; k++)
          for(s=0; s<nSpecies; s++)
          {
            p=particle[i][jend][k].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            } 
          }
      MPI_Send(&numP,1, MPI_INT,D->nextYrank,myrank, MPI_COMM_WORLD);    
    }    
    else if(rank%2==1) 
      MPI_Recv(&numP,1, MPI_INT,D->prevYrank,D->prevYrank, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    upP=(double *)malloc(totalData*sizeof(double ));             

    if(rank%2==0 && rank!=D->M-1)
    {    
      n=0;
      for(i=istart; i<iend; i++)
        for(k=kstart; k<kend; k++)
          for(s=0; s<nSpecies; s++)
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
              upP[n*sendData+9]=p->index;  
              upP[n*sendData+10]=(double)s;  
              upP[n*sendData+11]=(double)(p->core);  
              upP[n*sendData+12]=p->weight;  
              upP[n*sendData+13]=p->charge;  
              upP[n*sendData+14]=p->p1Old1;  
              upP[n*sendData+15]=p->p2Old1;  
              upP[n*sendData+16]=p->p3Old1;  
              upP[n*sendData+17]=p->oldX2;  
              upP[n*sendData+18]=p->oldY2-nySub;  
              upP[n*sendData+19]=p->oldZ2;  
              p=p->next;
              n++;
            }
          }
      MPI_Send(upP,totalData, MPI_DOUBLE,D->nextYrank, myrank, MPI_COMM_WORLD); 
    }

    else if(rank%2==1) 
    {    
      MPI_Recv(upP,totalData, MPI_DOUBLE,D->prevYrank,D->prevYrank, MPI_COMM_WORLD,&status);      
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
        New->index=upP[n*sendData+9];
        New->core=(int)(upP[n*sendData+11]);
        New->weight=upP[n*sendData+12];
        New->charge=upP[n*sendData+13];
        New->p1Old1=upP[n*sendData+14];
        New->p2Old1=upP[n*sendData+15];
        New->p3Old1=upP[n*sendData+16];
        New->oldX2=upP[n*sendData+17];
        New->oldY2=upP[n*sendData+18];
        New->oldZ2=upP[n*sendData+19];
      }
    }
    free(upP);
    MPI_Barrier(MPI_COMM_WORLD);

    //Odd -> evem
    if(rank%2==1 && rank!=D->M-1)
    {
      numP=0;
      for(i=istart; i<iend; i++)
        for(k=kstart; k<kend; k++)
          for(s=0; s<nSpecies; s++)
          {
            p=particle[i][jend][k].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            }
          } 
      MPI_Send(&numP,1, MPI_INT,D->nextYrank, myrank, MPI_COMM_WORLD);    
    }
    else if(rank%2==0 && rank!=0) 
      MPI_Recv(&numP,1, MPI_INT,D->prevYrank,D->prevYrank, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    upP=(double *)malloc(totalData*sizeof(double ));             

    if(rank%2==1 && rank!=D->M-1)
    {
      n=0;
      for(i=istart; i<iend; i++)
        for(k=kstart; k<kend; k++)
          for(s=0; s<nSpecies; s++)
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
              upP[n*sendData+9]=p->index;   
              upP[n*sendData+10]=(double)s;   
              upP[n*sendData+11]=(double)(p->core);   
              upP[n*sendData+12]=p->weight;  
              upP[n*sendData+13]=p->charge;  
              upP[n*sendData+14]=p->p1Old1;  
              upP[n*sendData+15]=p->p2Old1;  
              upP[n*sendData+16]=p->p3Old1;  
              upP[n*sendData+17]=p->oldX2;  
              upP[n*sendData+18]=p->oldY2-nySub;  
              upP[n*sendData+19]=p->oldZ2;  
              p=p->next;
              n++;
            }
        }
      MPI_Send(upP,totalData, MPI_DOUBLE,D->nextYrank, myrank, MPI_COMM_WORLD); 
    }
 
    else if(rank%2==0 && rank!=0) 
    {    
      MPI_Recv(upP,totalData, MPI_DOUBLE,D->prevYrank,D->prevYrank, MPI_COMM_WORLD,&status);      
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
        New->index=upP[n*sendData+9];
        New->core=(int)(upP[n*sendData+11]);
        New->weight=upP[n*sendData+12];
        New->charge=upP[n*sendData+13];
        New->p1Old1=upP[n*sendData+14];
        New->p2Old1=upP[n*sendData+15];
        New->p3Old1=upP[n*sendData+16];
        New->oldX2=upP[n*sendData+17];
        New->oldY2=upP[n*sendData+18];
        New->oldZ2=upP[n*sendData+19];
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(upP);

}


void MPI_TransferP_Yminus(Domain *D,Particle ***particle,int nSpecies)
{
   int i,j,k,n,s,numP,cnt,totalData,sendData=20,nxSub,nySub,nzSub;
  int istart,iend,jstart,jend,kstart,kend;
   int myrank, nTasks, rank;    
   double *btP;
   ptclList *p,*tmp,*New;
   MPI_Status status;         

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  nxSub=D->nxSub;       nySub=0;            nzSub=0;    
  istart=D->istart;     iend=D->iend;
  jstart=0;             jend=1;
  kstart=0;             kend=1;  
	if(D->dimension>1) {	jstart=D->jstart;   jend=D->jend;   nySub=D->nySub; } else ;
	if(D->dimension>2) {  kstart=D->kstart-1; kend=D->kend+1; nzSub=D->nzSub; } else ;

   rank=(myrank%(D->M*D->N))%D->M;

    //Even -> odd
    if(rank%2==0 && rank!=0)
    {
      numP=0;
      for(i=istart; i<iend; i++)
        for(k=kstart; k<kend; k++)
          for(s=0; s<nSpecies; s++)
          { 
            p=particle[i][jstart-1][k].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            } 
          }
      MPI_Send(&numP,1, MPI_INT,D->prevYrank, myrank, MPI_COMM_WORLD);    
    }    
    else if(rank%2==1 && rank!=D->M-1) 
      MPI_Recv(&numP,1, MPI_INT,D->nextYrank,D->nextYrank,MPI_COMM_WORLD,&status);        
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    btP=(double *)malloc(totalData*sizeof(double ));             
 
    if(rank%2==0 && rank!=0)
    {
      n=0;
      for(i=istart; i<iend; i++)
        for(k=kstart; k<kend; k++)
          for(s=0; s<nSpecies; s++)
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
              btP[n*sendData+9]=p->index;  
              btP[n*sendData+10]=(double)s;  
              btP[n*sendData+11]=(double)(p->core);  
              btP[n*sendData+12]=p->weight;  
              btP[n*sendData+13]=p->charge;  
              btP[n*sendData+14]=p->p1Old1;  
              btP[n*sendData+15]=p->p2Old1;  
              btP[n*sendData+16]=p->p3Old1;  
              btP[n*sendData+17]=p->oldX2;  
              btP[n*sendData+18]=p->oldY2;  
              btP[n*sendData+19]=p->oldZ2;  
              p=p->next;
              n++;
            }
          }
      MPI_Send(btP,totalData, MPI_DOUBLE,D->prevYrank,myrank, MPI_COMM_WORLD); 
    }
    else if(rank%2==1 && rank!=D->M-1) 
    {    
      MPI_Recv(btP,totalData, MPI_DOUBLE,D->nextYrank,D->nextYrank,MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
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
        New->index=btP[n*sendData+9];
        New->core=(int)(btP[n*sendData+11]);
        New->weight=btP[n*sendData+12];
        New->charge=btP[n*sendData+13];
        New->p1Old1=btP[n*sendData+14];
        New->p2Old1=btP[n*sendData+15];
        New->p3Old1=btP[n*sendData+16];
        New->oldX2=btP[n*sendData+17];
        New->oldY2=btP[n*sendData+18]+nySub;
        New->oldZ2=btP[n*sendData+19];
      }
    }
    free(btP);
    MPI_Barrier(MPI_COMM_WORLD);

    //Odd -> evem
    if(rank%2==1)
    {
      numP=0;
      for(i=istart; i<iend; i++)
        for(k=kstart; k<kend; k++)
          for(s=0; s<nSpecies; s++)
          {
            p=particle[i][jstart-1][k].head[s]->pt;
            while(p)   {
              p=p->next;
              numP++;
            }
          } 
      MPI_Send(&numP,1, MPI_INT,D->prevYrank, myrank, MPI_COMM_WORLD);    
    }
    else if(rank%2==0 && rank!=D->M-1) 
      MPI_Recv(&numP,1, MPI_INT,D->nextYrank,D->nextYrank,MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    btP=(double *)malloc(totalData*sizeof(double ));             

    if(rank%2==1)
    {
      n=0;
      for(i=istart; i<iend; i++)
        for(k=kstart; k<kend; k++)
          for(s=0; s<nSpecies; s++)
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
              btP[n*sendData+9]=p->index;  
              btP[n*sendData+10]=(double)s;  
              btP[n*sendData+11]=(double)(p->core);  
              btP[n*sendData+12]=p->weight;  
              btP[n*sendData+13]=p->charge;  
              btP[n*sendData+14]=p->p1Old1;  
              btP[n*sendData+15]=p->p2Old1;  
              btP[n*sendData+16]=p->p3Old1;  
              btP[n*sendData+17]=p->oldX2;  
              btP[n*sendData+18]=p->oldY2;  
              btP[n*sendData+19]=p->oldZ2;  
              p=p->next;
              n++;
            }
          }
      MPI_Send(btP,totalData, MPI_DOUBLE,D->prevYrank, myrank, MPI_COMM_WORLD); 
    }

    else if(rank%2==0 && rank!=D->M-1) 
    {    
      MPI_Recv(btP,totalData, MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
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
        New->index=btP[n*sendData+9];
        New->core=(int)(btP[n*sendData+11]);
        New->weight=btP[n*sendData+12];
        New->charge=btP[n*sendData+13];
        New->p1Old1=btP[n*sendData+14];
        New->p2Old1=btP[n*sendData+15];
        New->p3Old1=btP[n*sendData+16];
        New->oldX2=btP[n*sendData+17];
        New->oldY2=btP[n*sendData+18]+nySub;
        New->oldZ2=btP[n*sendData+19];
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(btP);

}

