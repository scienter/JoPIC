#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>


void boundary(Domain *D,External *Ext)
{
     FILE *out;
     int i,j,k,s,rank,rankX,rankY,rankZ,trackStart;
     int remainX,remainY,remainZ,subX,subY,subZ,tmpX,tmpY,tmpZ;
     int nx,ny,nz,nxSub,nySub,numdataUp,numdataBt;
	  long long numberData;
     int startj,startk,nxSub1D,nySub2D,nzSub3D;
     int minX,maxX,minY,maxY,minZ,maxZ;
     int myrank, nTasks,a;
     double ***memoryAsign(int nx, int ny, int nz);
     double ***memoryAsignJ(int nx, int ny, int nz);
     MPI_Status status;

     MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

     if(D->nx%D->resolX!=0 || D->ny%D->resolY!=0 || D->nz%D->resolZ!=0)  {
       printf("check nx, ny, nz. The remain of resol value is not zero.\n");
       exit(0);
     }  else ;
     nx=D->nx/D->resolX;
     ny=D->ny/D->resolY;
     nz=D->nz/D->resolZ;

     D->nxSub=nx/D->L;
     subX=D->nxSub;
     remainX=nx%D->L;
     minX=maxX=0;

     D->nySub=ny/D->M;
     subY=D->nySub;
     remainY=ny%D->M;
     minY=maxY=0;

     D->nzSub=nz/D->N;
     subZ=D->nzSub;
     remainZ=nz%D->N;
     minZ=maxZ=0;

     minX=maxX=0;
     for(rankX=0; rankX<D->L; rankX++)
     {
       if(rankX<remainX)   tmpX=subX+1;
       else                tmpX=subX;
       minX=maxX;
       maxX=minX+tmpX;

       minZ=maxZ=0;
       for(rankZ=0; rankZ<D->N; rankZ++)
       {
         if(rankZ<remainZ)   tmpZ=subZ+1;
         else                tmpZ=subZ;
         minZ=maxZ;
         maxZ=minZ+tmpZ;

         minY=maxY=0;
         for(rankY=0; rankY<D->M; rankY++)
         {
           if(rankY<remainY)   tmpY=subY+1;
           else                tmpY=subY;
           minY=maxY;
           maxY=minY+tmpY;

           rank=rankY+rankZ*D->M+rankX*(D->M*D->N);
           if(myrank==rank)
           {
              D->minXSub=minX*D->resolX;
              D->maxXSub=maxX*D->resolX;
              D->nxSub=tmpX*D->resolX;
              D->minYSub=minY*D->resolY;
              D->maxYSub=maxY*D->resolY;
              D->nySub=tmpY*D->resolY;
              D->minZSub=minZ*D->resolZ;
              D->maxZSub=maxZ*D->resolZ;
              D->nzSub=tmpZ*D->resolZ;
           }
         }
       }
     }
     D->minXSub+=D->minXDomain;    
     D->minYSub+=D->minYDomain;    
     D->minZSub+=D->minZDomain;    
     D->maxXSub+=D->minXDomain;    
     D->maxYSub+=D->minYDomain;    
     D->maxZSub+=D->minZDomain;    

     for(rank=0; rank<nTasks; rank++)
     {
       if(myrank==rank)
       {
//         if(D->dimension==1)
//         {
//           printf("rank=%d, minXSub=%d,maxSub=%d\n",myrank,D->minXSub,D->maxXSub);
//         }
//         else if(D->dimension==2)
//         {
//           printf("rank=%d, minXSub=%d,maxSub=%d,minYSub=%d,maxYSub=%d\n",myrank,D->minXSub,D->maxXSub,D->minYSub,D->maxYSub);
//         }
//         else if(D->dimension==3)
//         {
           printf("rank=%d, minXSub=%d,maxSub=%d,minYSub=%d,maxYSub=%d,minZSub=%d,maxZSub=%d\n",myrank,D->minXSub,D->maxXSub,D->minYSub,D->maxYSub,D->minZSub,D->maxZSub);
//         }
       }
     }

     //defining next and prev domain
     rankX=myrank/(D->M*D->N);
     rankZ=(myrank%(D->M*D->N))/D->M;
     rankY=(myrank%(D->M*D->N))%D->M;

     D->nextXrank=rankY+rankZ*D->M+(rankX+1)*(D->M*D->N);
     D->prevXrank=rankY+rankZ*D->M+(rankX-1)*(D->M*D->N);
     if(rankX==D->L-1 && rankX==0) {
       D->prevXrank=myrank;
       D->nextXrank=myrank;
     } else if(rankX==D->L-1) 
       D->nextXrank=rankY+rankZ*D->M;
     else if(rankX==0) 
       D->prevXrank=rankY+rankZ*D->M+(D->L-1)*(D->M*D->N);

     D->nextYrank=(rankY+1)+rankZ*D->M+rankX*(D->M*D->N);
     D->prevYrank=(rankY-1)+rankZ*D->M+rankX*(D->M*D->N);
     if(rankY==D->M-1 && rankY==0) {
       D->prevYrank=myrank;
       D->nextYrank=myrank;
     } else if(rankY==D->M-1) 
       D->nextYrank=rankZ*D->M+rankX*(D->M*D->N);
     else if(rankY==0) 
       D->prevYrank=(D->M-1)+rankZ*D->M+rankX*(D->M*D->N);


     D->nextZrank=rankY+(rankZ+1)*D->M+rankX*(D->M*D->N);
     D->prevZrank=rankY+(rankZ-1)*D->M+rankX*(D->M*D->N);
     if(rankZ==D->N-1 && rankZ==0) {
       D->prevZrank=myrank;
       D->nextZrank=myrank;
     } else if(rankZ==D->N-1) 
       D->nextZrank=rankY+rankX*(D->M*D->N);
     else if(rankZ==0) 
       D->prevZrank=rankY+(D->N-1)*D->M+rankX*(D->M*D->N);


     D->istart=2;
     D->iend=D->nxSub+2;
     D->jstart=0;
     D->jend=1;
     D->kstart=0;
     D->kend=1;
     if(D->dimension>1)  {
        D->jstart=2;
        D->jend=D->nySub+2;
     }
     if(D->dimension>2)  {
        D->kstart=2;
        D->kend=D->nzSub+2;
     }

     // Field setting
     nxSub1D=D->nxSub+5;
     nySub2D=1;
     nzSub3D=1;
     if(D->dimension>1)  
       nySub2D=D->nySub+5;
     if(D->dimension>2)  
       nzSub3D=D->nzSub+5;

	  D->F=memoryAsign(nxSub1D,nySub2D,nzSub3D);
	  D->Rho=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     D->RhoPair=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     D->RhoNoPair=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     D->Phi=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     D->PhiOld=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     D->Den=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     D->DenOld=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     D->CurX=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     D->CurY=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     D->CurZ=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     D->Ax=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     D->Ay=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     D->Az=memoryAsign(nxSub1D,nySub2D,nzSub3D);

     switch(D->fieldType)  {
     case Split :
       D->Ex=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Bx=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->ExNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->BxNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Pr=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Pl=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Sr=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Sl=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->PrC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->PlC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->SrC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->SlC=memoryAsign(nxSub1D,nySub2D,nzSub3D);		 
       D->Jx=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->Jy=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->Jz=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->Jz=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       break;

     case Yee :
     case Pukhov :
       D->Ex=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Ey=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Ez=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Bx=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->By=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Bz=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->BxNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->ByNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->BzNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Jx=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->Jy=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->Jz=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       break;

     default :
       printf("what field type?\n");
     }

     //PML
     D->frr=(double *)malloc(nzSub3D*sizeof(double ));
     D->bkr=(double *)malloc(nzSub3D*sizeof(double ));
     D->upr=(double *)malloc(nySub2D*sizeof(double ));
     D->dnr=(double *)malloc(nySub2D*sizeof(double ));
     D->rtr=(double *)malloc(nxSub1D*sizeof(double ));
     D->ltr=(double *)malloc(nxSub1D*sizeof(double ));
     D->frd=(double *)malloc(nzSub3D*sizeof(double ));
     D->bkd=(double *)malloc(nzSub3D*sizeof(double ));
     D->upd=(double *)malloc(nySub2D*sizeof(double ));
     D->dnd=(double *)malloc(nySub2D*sizeof(double ));
     D->rtd=(double *)malloc(nxSub1D*sizeof(double ));
     D->ltd=(double *)malloc(nxSub1D*sizeof(double ));

    //share
    D->shareX=3;  D->shareY=0;  D->shareZ=0;
    if(D->dimension>1)  D->shareY=3; else ;
    if(D->dimension>2)  D->shareZ=3; else ;
    D->shareF = (double ****)malloc(6*sizeof(double ***));



  // Particle setting
  nxSub1D=D->nxSub+3;
  nySub2D=1;
  nzSub3D=1;
  startj=0;
  startk=0;
  if(D->dimension>1) {
    nySub2D=D->nySub+3;
    startj=D->jstart-1;
  }
  if(D->dimension>2) { 
    nzSub3D=D->nzSub+3;
    startk=D->kstart-1;
  }

  D->index=1;
  D->particle = (Particle ***)malloc((nxSub1D)*sizeof(Particle **));
  for(i=0; i<nxSub1D; i++) {
    D->particle[i] = (Particle **)malloc((nySub2D)*sizeof(Particle *));
    for(j=0; j<nySub2D; j++) 
      D->particle[i][j] = (Particle *)malloc((nzSub3D)*sizeof(Particle ));
  }

  // setting up particle's pointer
  if(D->dimension==1)     {
    j=k=0;
    //i starts at 0 because of boost frame
    for(i=0; i<D->iend+1; i++)     {
      D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
      for(s=0; s<D->nSpecies; s++)      {
        D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
        D->particle[i][j][k].head[s]->pt = NULL;
      }
    }
  } else if(D->dimension==2)   {
    k=0;
    for(i=0; i<D->iend+1; i++)	//i starts at 0 because of boost frame
      for(j=D->jstart-1; j<D->jend+1; j++)   {
        D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
        for(s=0; s<D->nSpecies; s++)           {
          D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
          D->particle[i][j][k].head[s]->pt = NULL;
        }
      }
  } else if(D->dimension==3)     {
    for(i=0; i<D->iend+1; i++)	//i starts at 0 because of boost frame
      for(j=D->jstart-1; j<D->jend+1; j++) {
        for(k=D->kstart-1; k<D->kend+1; k++)           {
          D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
          for(s=0; s<D->nSpecies; s++)             {
            D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
            D->particle[i][j][k].head[s]->pt = NULL;
          }
        }
      }
  }

  //Making track Particle memory
	D->track = (Particle ***)malloc((nxSub1D)*sizeof(Particle **));
	for(i=0; i<nxSub1D; i++) {
		D->track[i] = (Particle **)malloc((nySub2D)*sizeof(Particle *));
		for(j=0; j<nySub2D; j++)
			D->track[i][j] = (Particle *)malloc((nzSub3D)*sizeof(Particle ));
	}
	if(D->dimension==1)     {
		j=k=0;
		for(i=0; i<D->iend+1; i++)     {
			D->track[i][j][k].head = (ptclHead **)malloc(1*sizeof(ptclHead *));
			D->track[i][j][k].head[0] = (ptclHead *)malloc(sizeof(ptclHead));
			D->track[i][j][k].head[0]->pt = NULL;
		}
	} else if(D->dimension==2)   {
  	k=0;
    for(i=0; i<D->iend+1; i++) //i starts at 0 because of boost frame
		  for(j=D->jstart-1; j<D->jend+1; j++)   {
			  D->track[i][j][k].head = (ptclHead **)malloc(1*sizeof(ptclHead *));
			  D->track[i][j][k].head[0] = (ptclHead *)malloc(sizeof(ptclHead));
			  D->track[i][j][k].head[0]->pt = NULL;
		  }
	} else if(D->dimension==3)     {
		for(i=0; i<D->iend+1; i++) //i starts at 0 because of boost frame
			for(j=D->jstart-1; j<D->jend+1; j++) {
				for(k=D->kstart-1; k<D->kend+1; k++)           {
					D->track[i][j][k].head = (ptclHead **)malloc(1*sizeof(ptclHead *));
					D->track[i][j][k].head[0] = (ptclHead *)malloc(sizeof(ptclHead));
					D->track[i][j][k].head[0]->pt = NULL;
				}
			}
	}

	int dataNum;
	D->trackStep=100;
	dataNum=D->trackStep*D->idNums*6;
	D->trackParticle=(double *)malloc(dataNum*sizeof(double ));
	for(i=0; i<dataNum; i++) {
		D->trackParticle[i]=0.0;
	}


}

double ***memoryAsign(int nx, int ny, int nz)
{
   int i,j,k;
   double ***field;

   field = (double ***)malloc((nx)*sizeof(double **));
   for(i=0; i<nx; i++)   {
     field[i] = (double **)malloc((ny)*sizeof(double *));
     for(j=0; j<ny; j++)
       field[i][j] = (double *)malloc((nz)*sizeof(double ));
   }
   
   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++)
         field[i][j][k]=0.0;

   return field;
}

double ***memoryAsignJ(int nx, int ny, int nz)
{
   int i,j,k;
   double ***current;

   current = (double ***)malloc((nx)*sizeof(double **));
   for(i=0; i<nx; i++)   {
     current[i] = (double **)malloc((ny)*sizeof(double *));
     for(j=0; j<ny; j++)
       current[i][j] = (double *)malloc((nz)*sizeof(double ));
   }
   
   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++)
         current[i][j][k]=0.0;

   return current;
}
