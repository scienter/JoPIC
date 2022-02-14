#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>


void MPI_TransferP_Period_X(Domain *D,int istart,int iend,int jstart,int jend)
{
  int i,j,n,s,numP,cnt,totalData,sendData=15,nxSub,nySub;
  int myrank, nTasks, rank;
  Particle **particle;
  particle=D->particle;     
  double *data;
  ptclList *p,*New;

  MPI_Status status;         
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  nxSub=D->nxSub; nySub=D->nySub;
  rank=myrank/D->M;

  //Right to Left
	numP=0;
  for(j=jstart; j<jend; j++)
   	for(s=0; s<D->nSpecies; s++) {
      	p=particle[iend][j].head[s]->pt;
        while(p)   {
         	p=p->next;
          numP++;
       	} 
  		}
	if(nTasks>1) {
		if(myrank==nTasks-1) 
      	MPI_Send(&numP,1, MPI_INT,0, myrank, MPI_COMM_WORLD); 
   	else if(myrank==0) 
      	MPI_Recv(&numP,1,MPI_INT,nTasks-1,nTasks-1,MPI_COMM_WORLD,&status);      
    	MPI_Barrier(MPI_COMM_WORLD);    
	} else ;

	totalData=numP*sendData;
	data=(double *)malloc(totalData*sizeof(double ));             

	if(nTasks==1) {
		n=0;
      for(j=jstart; j<jend; j++)
        for(s=0; s<D->nSpecies; s++) {
          p=particle[iend][j].head[s]->pt;
          while(p) {
            data[n*sendData+0]=iend+p->z-nxSub;     
            data[n*sendData+1]=p->oldZ-nxSub;  
            data[n*sendData+2]=p->x;  
            data[n*sendData+3]=p->oldX;  
            data[n*sendData+4]=p->y;  
            data[n*sendData+5]=p->oldY;  
            data[n*sendData+6]=p->pz;  
            data[n*sendData+7]=p->px;  
            data[n*sendData+8]=p->py;  
            data[n*sendData+9]=(double)(p->index);  
            data[n*sendData+10]=(double)s;  
            data[n*sendData+11]=(double)(p->core);  
            data[n*sendData+12]=p->weight;  
            data[n*sendData+13]=p->charge;  
            data[n*sendData+14]=(double)j;  
            p=p->next;
            n++;
          }
        }
    for(n=0; n<numP; n++) {
    	i=(int)(data[n*sendData+0]);
			j=(int)(data[n*sendData+14]);
			s=(int)(data[n*sendData+10]);
	    New = (ptclList *)malloc(sizeof(ptclList)); 
    	New->next = particle[i][j].head[s]->pt;
	    particle[i][j].head[s]->pt = New;             
    	New->z   =data[n*sendData+0]-i;     
			New->oldZ=data[n*sendData+1];     
   		New->x   =data[n*sendData+2];     
			New->oldX=data[n*sendData+3];     
			New->y   =data[n*sendData+4];     
			New->oldY=data[n*sendData+5];     
      New->pz  =data[n*sendData+6];
      New->px  =data[n*sendData+7];
      New->py  =data[n*sendData+8];
      New->index=(int)(data[n*sendData+9]);
      New->core =(int)(data[n*sendData+11]);
      New->weight=data[n*sendData+12];
      New->charge=data[n*sendData+13];
    }
	} else {
		if(myrank==nTasks-1) {
      n=0;
	    for(j=jstart; j<jend; j++)
   		  for(s=0; s<D->nSpecies; s++) {
        	p=particle[iend][j].head[s]->pt;
        	while(p) {
          	data[n*sendData+0]=iend+p->z-nxSub;     
          	data[n*sendData+1]=p->oldZ-nxSub;  
          	data[n*sendData+2]=p->x;  
          	data[n*sendData+3]=p->oldX;  
          	data[n*sendData+4]=p->y;  
          	data[n*sendData+5]=p->oldY;  
          	data[n*sendData+6]=p->pz;  
		        data[n*sendData+7]=p->px;  
      	    data[n*sendData+8]=p->py;  
		        data[n*sendData+9]=(double)(p->index);  
      	    data[n*sendData+10]=(double)s;  
		        data[n*sendData+11]=(double)(p->core);  
      	    data[n*sendData+12]=p->weight;  
		        data[n*sendData+13]=p->charge;  
      	    data[n*sendData+14]=(double)j;  
          	p=p->next;
		        n++;
      		}
      	}
			MPI_Send(data,totalData, MPI_DOUBLE,D->nextXrank, myrank, MPI_COMM_WORLD); 
    }
		else if(myrank==0) {    
      MPI_Recv(data,totalData, MPI_DOUBLE,nTasks-1,nTasks-1, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++) {
      	i=(int)(data[n*sendData+0]);
			  j=(int)(data[n*sendData+14]);
			  s=(int)(data[n*sendData+10]);
		    New = (ptclList *)malloc(sizeof(ptclList)); 
      	New->next = particle[i][j].head[s]->pt;
		    particle[i][j].head[s]->pt = New;             
      	New->z   =data[n*sendData+0]-i;     
	      New->oldZ=data[n*sendData+1];     
    		New->x   =data[n*sendData+2];     
			  New->oldX=data[n*sendData+3];     
			  New->y   =data[n*sendData+4];     
			  New->oldY=data[n*sendData+5];     
      	New->pz  =data[n*sendData+6];
      	New->px  =data[n*sendData+7];
      	New->py  =data[n*sendData+8];
      	New->index=(int)(data[n*sendData+9]);
      	New->core =(int)(data[n*sendData+11]);
      	New->weight=data[n*sendData+12];
      	New->charge=data[n*sendData+13];
      }
    }   
    MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);

  // Left -> Right
	numP=0;
  for(i=0; i<istart; i++)
   	for(j=jstart; j<jend; j++)
      for(s=0; s<D->nSpecies; s++)  { 
       	p=particle[i][j].head[s]->pt;
        while(p)   {
          p=p->next;
          numP++;
        } 
      }
	if(nTasks>1) {
		if(myrank==0) 
      MPI_Send(&numP,1, MPI_INT,nTasks-1, myrank, MPI_COMM_WORLD);
    else if(myrank==nTasks-1) 
      MPI_Recv(&numP,1, MPI_INT,0,0,MPI_COMM_WORLD,&status);        
		MPI_Barrier(MPI_COMM_WORLD);    
	} else ;

	totalData=numP*sendData;
  data=(double *)malloc(totalData*sizeof(double ));             
 
	if(nTasks==1) {
    n=0;
    for(i=0; i<istart; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<D->nSpecies; s++) {
          p=particle[i][j].head[s]->pt;
          while(p)   {
            data[n*sendData+0]=p->z+i;     
            data[n*sendData+1]=p->oldZ+nxSub;  
            data[n*sendData+2]=p->x;  
            data[n*sendData+3]=p->oldX;  
            data[n*sendData+4]=p->y;  
            data[n*sendData+5]=p->oldY;  
            data[n*sendData+6]=p->pz;  
            data[n*sendData+7]=p->px;  
            data[n*sendData+8]=p->py;  
            data[n*sendData+9]=(double)(p->index);  
            data[n*sendData+10]=(double)s;  
            data[n*sendData+11]=(double)(p->core);  
            data[n*sendData+12]=p->weight;  
            data[n*sendData+13]=p->charge;  
            data[n*sendData+14]=(double)j;  
            p=p->next;
            n++;
          }
        }
    for(n=0; n<numP; n++) {
      i=(int)(data[n*sendData+0]);
      j=(int)(data[n*sendData+14]);
      s=(int)(data[n*sendData+10]);

      New = (ptclList *)malloc(sizeof(ptclList)); 
      New->next = particle[i+nxSub][j].head[s]->pt;
      particle[i+nxSub][j].head[s]->pt = New;             
      New->z   =data[n*sendData+0]-i;     
      New->oldZ=data[n*sendData+1];     
      New->x   =data[n*sendData+2];     
      New->oldX=data[n*sendData+3];     
      New->y   =data[n*sendData+4];     
      New->oldY=data[n*sendData+5];     
      New->pz  =data[n*sendData+6];
      New->px  =data[n*sendData+7];
      New->py  =data[n*sendData+8];
      New->index=(int)(data[n*sendData+9]);
      New->core =(int)(data[n*sendData+11]);
      New->weight=data[n*sendData+12];
      New->charge=data[n*sendData+13];
    }
	} else {
		if(myrank==0) {
      n=0;
			  for(i=0; i<istart; i++)
        	for(j=jstart; j<jend; j++)
        		for(s=0; s<D->nSpecies; s++) {
          		p=particle[i][j].head[s]->pt;
          		while(p)   {
            		data[n*sendData+0]=p->z+i;     
			          data[n*sendData+1]=p->oldZ;  
        		    data[n*sendData+2]=p->x;  
		 			      data[n*sendData+3]=p->oldX;  
		          	data[n*sendData+4]=p->y;  
      	      	data[n*sendData+5]=p->oldY;  
          			data[n*sendData+6]=p->pz;  
		          	data[n*sendData+7]=p->px;  
      	      	data[n*sendData+8]=p->py;  
          			data[n*sendData+9]=(double)(p->index);  
		          	data[n*sendData+10]=(double)s;  
      	      	data[n*sendData+11]=(double)(p->core);  
		          	data[n*sendData+12]=p->weight;  
      	      	data[n*sendData+13]=p->charge;  
		          	data[n*sendData+14]=(double)j;  
      	      	p=p->next;
          			n++;
		          }
					}
			MPI_Send(data,totalData, MPI_DOUBLE,nTasks-1, myrank, MPI_COMM_WORLD); 
		}
		else if(myrank==nTasks-1) {
      MPI_Recv(data,totalData, MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);      
	    for(n=0; n<numP; n++) {
   	    i=(int)(data[n*sendData+0]);
        j=(int)(data[n*sendData+14]);
	      s=(int)(data[n*sendData+10]);
	
   		  New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i+nxSub][j].head[s]->pt;
	      particle[i+nxSub][j].head[s]->pt = New;             
   	   	New->z   =data[n*sendData+0]-i;     
        New->oldZ=data[n*sendData+1]+nxSub;     
	    	New->x   =data[n*sendData+2];     
  			New->oldX=data[n*sendData+3];     
    		New->y   =data[n*sendData+4];     
	   	  New->oldY=data[n*sendData+5];     
   		  New->pz  =data[n*sendData+6];
    		New->px  =data[n*sendData+7];
	    	New->py  =data[n*sendData+8];
   		  New->index=(int)(data[n*sendData+9]);
	  		New->core =(int)(data[n*sendData+11]);
  			New->weight=data[n*sendData+12];
    		New->charge=data[n*sendData+13];
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);
}

void particleShareX(Domain D)
{
  int istart,iend,jstart,jend,kstart,kend;
  void MPI_TransferP_Xplus();
  void MPI_TransferP_Xminus();
 
  Particle **particle;
  particle=D.particle; 

  istart=D.istart;    iend=D.iend;
  jstart=D.jstart;    jend=D.jend;

  MPI_TransferP_Xplus(&D,istart,iend,jstart-1,jend+1);
  MPI_TransferP_Xminus(&D,istart,iend,jstart-1,jend+1);
}

void MPI_TransferP_Xplus(Domain *D,int istart,int iend,int jstart,int jend)
{
  int i,j,n,s,numP,cnt,totalData,sendData=15,nxSub,nySub;
  int myrank, nTasks, rank;
  Particle **particle;
  particle=D->particle;     
  double *upP;
  ptclList *p,*New;
  MPI_Status status;         

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  nxSub=D->nxSub; nySub=D->nySub;

  rank=myrank/D->M;

  //Even -> odd
  if(rank%2==0 && rank!=D->L-1) {
    numP=0;
    for(j=jstart; j<jend; j++)
      for(s=0; s<D->nSpecies; s++) {
        p=particle[iend][j].head[s]->pt;
        while(p)   {
          p=p->next;
          numP++;
        } 
      }
    MPI_Send(&numP,1, MPI_INT,D->nextXrank, myrank, MPI_COMM_WORLD); 
  }
  else if(rank%2==1) 
    MPI_Recv(&numP,1,MPI_INT,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);      
  MPI_Barrier(MPI_COMM_WORLD);    

  totalData=numP*sendData;
  upP=(double *)malloc(totalData*sizeof(double ));             

  if(rank%2==0 && rank!=D->L-1) {    
    n=0;
    for(j=jstart; j<jend; j++)
      for(s=0; s<D->nSpecies; s++) {
        p=particle[iend][j].head[s]->pt;
        while(p) {
          upP[n*sendData+0]=iend+p->z-nxSub;     
          upP[n*sendData+1]=p->oldZ-nxSub;  
          upP[n*sendData+2]=p->x;  
          upP[n*sendData+3]=p->oldX;  
          upP[n*sendData+4]=p->y;  
          upP[n*sendData+5]=p->oldY;  
          upP[n*sendData+6]=p->pz;  
          upP[n*sendData+7]=p->px;  
          upP[n*sendData+8]=p->py;  
          upP[n*sendData+9]=(double)(p->index);  
          upP[n*sendData+10]=(double)s;  
          upP[n*sendData+11]=(double)(p->core);  
          upP[n*sendData+12]=p->weight;  
          upP[n*sendData+13]=p->charge;  
          upP[n*sendData+14]=(double)j;  
          p=p->next;
          n++;
        }
      }
    MPI_Send(upP,totalData, MPI_DOUBLE,D->nextXrank, myrank, MPI_COMM_WORLD); 
  }

  else if(rank%2==1) {    
    MPI_Recv(upP,totalData, MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);      
    for(n=0; n<numP; n++) {
      i=(int)(upP[n*sendData+0]);
      j=(int)(upP[n*sendData+14]);
      s=(int)(upP[n*sendData+10]);
      New = (ptclList *)malloc(sizeof(ptclList)); 
      New->next = particle[i][j].head[s]->pt;
      particle[i][j].head[s]->pt = New;             
      New->z=upP[n*sendData+0]-i;     
      New->oldZ=upP[n*sendData+1];     
      New->x=upP[n*sendData+2];     
      New->oldX=upP[n*sendData+3];     
      New->y=upP[n*sendData+4];     
      New->oldY=upP[n*sendData+5];     
      New->pz=upP[n*sendData+6];
      New->px=upP[n*sendData+7];
      New->py=upP[n*sendData+8];
      New->index=(int)(upP[n*sendData+9]);
      New->core=(int)(upP[n*sendData+11]);
      New->weight=upP[n*sendData+12];
      New->charge=upP[n*sendData+13];
    }
  }   
  free(upP);
  MPI_Barrier(MPI_COMM_WORLD);

  //Odd -> evem
  if(rank%2==1 && rank!=D->L-1) {
    numP=0;
    for(j=jstart; j<jend; j++)
      for(s=0; s<D->nSpecies; s++) {
        p=particle[iend][j].head[s]->pt;
        while(p)   {
          p=p->next;
          numP++;
        }
      } 
    MPI_Send(&numP,1, MPI_INT,D->nextXrank, myrank, MPI_COMM_WORLD);
  }
  else if(rank%2==0 && rank!=0) 
    MPI_Recv(&numP,1, MPI_INT,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);      
  MPI_Barrier(MPI_COMM_WORLD);    

  totalData=numP*sendData;
  upP=(double *)malloc(totalData*sizeof(double ));             

  if(rank%2==1 && rank!=D->L-1) {
    n=0;
    for(j=jstart; j<jend; j++)
      for(s=0; s<D->nSpecies; s++)  {
        p=particle[iend][j].head[s]->pt;
        while(p)   {
          upP[n*sendData+0]=p->z+iend-nxSub;     
          upP[n*sendData+1]=p->oldZ-nxSub;  
          upP[n*sendData+2]=p->x;  
          upP[n*sendData+3]=p->oldX;  
          upP[n*sendData+4]=p->y;     
          upP[n*sendData+5]=p->oldY;  
          upP[n*sendData+6]=p->pz;  
          upP[n*sendData+7]=p->px;  
          upP[n*sendData+8]=p->py;  
          upP[n*sendData+9]=(double)(p->index);   
          upP[n*sendData+10]=(double)s;   
          upP[n*sendData+11]=(double)(p->core);   
          upP[n*sendData+12]=p->weight;  
          upP[n*sendData+13]=p->charge;  
          upP[n*sendData+14]=(double)j;  
          p=p->next;
          n++;
        }
      }
    MPI_Send(upP,totalData, MPI_DOUBLE,D->nextXrank, myrank, MPI_COMM_WORLD); 
  }
 
  else if(rank%2==0 && rank!=0) {    
    MPI_Recv(upP,totalData, MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);      
    for(n=0; n<numP; n++) {
      i=(int)(upP[n*sendData+0]);
      j=(int)(upP[n*sendData+14]);
      s=(int)(upP[n*sendData+10]);
      New = (ptclList *)malloc(sizeof(ptclList)); 
      New->next = particle[i][j].head[s]->pt;
      particle[i][j].head[s]->pt = New;             
      New->z=upP[n*sendData+0]-i;     
      New->oldZ=upP[n*sendData+1];     
      New->x=upP[n*sendData+2];     
      New->oldX=upP[n*sendData+3];     
      New->y=upP[n*sendData+4]; 
      New->oldY=upP[n*sendData+5];     
      New->pz=upP[n*sendData+6];
      New->px=upP[n*sendData+7];
      New->py=upP[n*sendData+8];
      New->index=(int)(upP[n*sendData+9]);
      New->core=(int)(upP[n*sendData+11]);
      New->weight=upP[n*sendData+12];
      New->charge=upP[n*sendData+13];
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  free(upP);
}

void MPI_TransferP_Xminus(Domain *D,int istart,int iend,int jstart,int jend)
{
  int i,j,n,s,numP,cnt,totalData,sendData=15,nxSub,nySub;
  int myrank, nTasks, rank;
  Particle **particle;
  particle=D->particle;     
  double *btP;
  ptclList *p,*New;

  MPI_Status status;         
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  nxSub=D->nxSub;   nySub=D->nySub;
  rank=myrank/D->M;

  //Even -> odd
  if(rank%2==0 && rank!=0) {
    numP=0;
    for(i=0; i<istart; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<D->nSpecies; s++)  { 
          p=particle[i][j].head[s]->pt;
          while(p)   {
            p=p->next;
            numP++;
          } 
        }
    MPI_Send(&numP,1, MPI_INT,D->prevXrank, myrank, MPI_COMM_WORLD);
  }    
  else if(rank%2==1 && rank!=D->L-1) 
    MPI_Recv(&numP,1, MPI_INT,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);        
  MPI_Barrier(MPI_COMM_WORLD);    

  totalData=numP*sendData;
  btP=(double *)malloc(totalData*sizeof(double ));             
 
  if(rank%2==0 && rank!=0) {
    n=0;
    for(i=0; i<istart; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<D->nSpecies; s++) {
          p=particle[i][j].head[s]->pt;
          while(p)   {
            btP[n*sendData+0]=p->z+i;     
            btP[n*sendData+1]=p->oldZ;  
            btP[n*sendData+2]=p->x;  
            btP[n*sendData+3]=p->oldX;  
            btP[n*sendData+4]=p->y;  
            btP[n*sendData+5]=p->oldY;  
            btP[n*sendData+6]=p->pz;  
            btP[n*sendData+7]=p->px;  
            btP[n*sendData+8]=p->py;  
            btP[n*sendData+9]=(double)(p->index);  
            btP[n*sendData+10]=(double)s;  
            btP[n*sendData+11]=(double)(p->core);  
            btP[n*sendData+12]=p->weight;  
            btP[n*sendData+13]=p->charge;  
            btP[n*sendData+14]=(double)j;  
            p=p->next;
            n++;
          }
        }
    MPI_Send(btP,totalData, MPI_DOUBLE,D->prevXrank, myrank, MPI_COMM_WORLD); 
  }
  else if(rank%2==1 && rank!=D->L-1) {    
    MPI_Recv(btP,totalData, MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);      
    for(n=0; n<numP; n++) {
      i=(int)(btP[n*sendData+0]);
      j=(int)(btP[n*sendData+14]);
      s=(int)(btP[n*sendData+10]);

      New = (ptclList *)malloc(sizeof(ptclList)); 
      New->next = particle[i+nxSub][j].head[s]->pt;
      particle[i+nxSub][j].head[s]->pt = New;             
      New->z=btP[n*sendData+0]-i;     
      New->oldZ=btP[n*sendData+1]+nxSub;     
      New->x=btP[n*sendData+2];     
      New->oldX=btP[n*sendData+3];     
      New->y=btP[n*sendData+4];     
      New->oldY=btP[n*sendData+5];     
      New->pz=btP[n*sendData+6];
      New->px=btP[n*sendData+7];
      New->py=btP[n*sendData+8];
      New->index=(int)(btP[n*sendData+9]);
      New->core=(int)(btP[n*sendData+11]);
      New->weight=btP[n*sendData+12];
      New->charge=btP[n*sendData+13];
    }
  }
  free(btP);
  MPI_Barrier(MPI_COMM_WORLD);

  //Odd -> evem	
  if(rank%2==1)  {
    numP=0;
    for(i=0; i<istart; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<D->nSpecies; s++) {
          p=particle[i][j].head[s]->pt;
          while(p)   {
            p=p->next;
            numP++;
          }
        } 
    MPI_Send(&numP,1, MPI_INT,D->prevXrank,myrank, MPI_COMM_WORLD);
  }
  else if(rank%2==0 && rank!=D->L-1) 
    MPI_Recv(&numP,1, MPI_INT,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);      
  
  MPI_Barrier(MPI_COMM_WORLD);    

  totalData=numP*sendData;
  btP=(double *)malloc(totalData*sizeof(double ));            
 
  if(rank%2==1) {
    n=0;
    for(i=0; i<istart; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<D->nSpecies; s++) {
          p=particle[i][j].head[s]->pt;
          while(p) {
            btP[n*sendData+0]=p->z+i;     
            btP[n*sendData+1]=p->oldZ;  
            btP[n*sendData+2]=p->x; 
            btP[n*sendData+3]=p->oldX;  
            btP[n*sendData+4]=p->y;  
            btP[n*sendData+5]=p->oldY;  
            btP[n*sendData+6]=p->pz;  
            btP[n*sendData+7]=p->px;  
            btP[n*sendData+8]=p->py;  
            btP[n*sendData+9]=(double)(p->index);  
            btP[n*sendData+10]=(double)s;  
            btP[n*sendData+11]=(double)(p->core);  
            btP[n*sendData+12]=p->weight;  
            btP[n*sendData+13]=p->charge;  
            btP[n*sendData+14]=(double)j;
            p=p->next;
            n++;
          }
        }
    MPI_Send(btP,totalData, MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD); 
  }

  else if(rank%2==0 && rank!=D->L-1)  {    
    MPI_Recv(btP,totalData, MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);     
 
    for(n=0; n<numP; n++) {
      i=(int)(btP[n*sendData+0]);
      j=(int)(btP[n*sendData+14]);
      s=(int)(btP[n*sendData+10]);
      New = (ptclList *)malloc(sizeof(ptclList)); 
      New->next = particle[i+nxSub][j].head[s]->pt;
      particle[i+nxSub][j].head[s]->pt = New;             
      New->z=btP[n*sendData+0]-i;     
      New->oldZ=btP[n*sendData+1]+nxSub;     
      New->x=btP[n*sendData+2];     
      New->oldX=btP[n*sendData+3];     
      New->y=btP[n*sendData+4];     
      New->oldY=btP[n*sendData+5];     
      New->pz=btP[n*sendData+6];
      New->px=btP[n*sendData+7];
      New->py=btP[n*sendData+8];
      New->index=(int)(btP[n*sendData+9]);
      New->core=(int)(btP[n*sendData+11]);
      New->weight=btP[n*sendData+12];
      New->charge=btP[n*sendData+13];
    }
  }
  free(btP);
  MPI_Barrier(MPI_COMM_WORLD);
}
