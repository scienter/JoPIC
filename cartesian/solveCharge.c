#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void solveF_Split(Domain *D,double ***Pr,double ***Pl,double ***Sr,double ***Sl,double ***Ex,double ***Bx,int half,int non_half);

void solveF(Domain *D,double ***Pr,double ***Pl,double ***Sr,double ***Sl,double ***Ex,double ***Bx,int half,int non_half)
{
   switch (D->fieldType) {
   case Yee :
   case Pukhov :
     //solveF_YeePukhov(&D);
     break;
   case Split :
     solveF_Split(D,D->Pr,D->Pl,D->Sr,D->Sl,D->ExNow,D->BxNow,half,non_half);
     break;
   }
}

void solveF_Split(Domain *D,double ***Pr,double ***Pl,double ***Sr,double ***Sl,double ***Ex,double ***Bx,int half,int non_half)
{
  	int i,j,k,s,rangeY,rangeZ;
  	int istart,iend,jstart,jend,kstart,kend;
  	double invDx,invDy,invDz;
  	LoadList *LL;

  	istart=D->istart;  iend=D->iend;
  	jstart=D->jstart;  jend=D->jend;
  	kstart=D->kstart;  kend=D->kend;

  	invDx=1.0/D->dx;
  	invDy=1.0/D->dy;
  	invDz=1.0/D->dz;

  	// initialize RhoNoPair
  	rangeY=1; rangeZ=1;
  	if(D->dimension>1)  rangeY=D->jend+3; else ;
  	if(D->dimension>2)  rangeZ=D->kend+3; else ;
		for(i=0; i<iend+3; i++)
			for(j=0; j<rangeY; j++)
				for(k=0; k<rangeZ; k++)
  	      		D->RhoNoPair[i][j][k]=0.0;  

  	// calculation density
	LL=D->loadList;      
  	s=0;
  	while(LL->next)  {
  	   solveCharge(D,LL,D->RhoNoPair,s,istart,iend,1.0,half,non_half);
  	  	LL=LL->next;
  	  	s++;
  	}
  	D->shareF[0]=D->RhoNoPair;
  	MPI_TransferJ_Xplus(D,1);
  	MPI_TransferJ_Xminus(D,1);
  	MPI_TransferJ_Yplus(D,1);
  	MPI_TransferJ_Yminus(D,1);
  	MPI_TransferJ_Zplus(D,1);
  	MPI_TransferJ_Zminus(D,1);
  	D->RhoNoPair=D->shareF[0];

  	// solve F
  	switch (D->dimension) {
  	case 2 :
    	k=0;
	  	for(i=istart; i<iend; i++)
	  		for(j=jstart; j<jend; j++)
        	{
          	D->F[i][j][k]=invDx*(Ex[i][j][k]-Ex[i-1][j][k]);
            	+0.5*invDy*(Pr[i][j][k]+Pr[i-1][j][k]-Pr[i][j-1][k]-Pr[i-1][j-1][k]+Pl[i][j][k]+Pl[i-1][j][k]-Pl[i][j-1][k]-Pl[i-1][j-1][k])
            	-2.0*M_PI*(D->RhoNoPair[i][j][k]+D->RhoPair[i][j][k]);
        	}
    	break;
  	case 3 :
	  	for(i=istart; i<iend; i++)
	  		for(j=jstart; j<jend; j++)
	  			for(k=jstart; k<kend; k++) 
        		{
          		D->F[i][j][k]=invDx*(Ex[i][j][k]-Ex[i-1][j][k])
            		+0.5*invDy*(Pr[i][j][k]+Pr[i-1][j][k]-Pr[i][j-1][k]-Pr[i-1][j-1][k]+Pl[i][j][k]+Pl[i-1][j][k]-Pl[i][j-1][k]-Pl[i-1][j-1][k])
            		+0.5*invDz*(Sr[i][j][k]+Sr[i-1][j][k]-Sr[i][j][k-1]-Sr[i-1][j][k-1]+Sl[i][j][k]+Sl[i-1][j][k]-Sl[i][j][k-1]-Sl[i-1][j][k-1])
            		-2.0*M_PI*(D->RhoNoPair[i][j][k]+D->RhoPair[i][j][k]);
        		}
    	break;
  	}

  	D->shareF[0]=D->F;
  	MPI_TransferF_Xminus(D,1);
  	MPI_TransferF_Xplus(D,1);
  	MPI_TransferF_Yminus(D,1);
  	MPI_TransferF_Yminus(D,1);
  	MPI_TransferF_Zplus(D,1);
  	MPI_TransferF_Zplus(D,1);
  	if(D->Period==ON)  {
  	  	MPI_TransferF_Period_X(D,1);
  	  	MPI_TransferF_Period_Y(D,1);
  	  	MPI_TransferF_Period_Z(D,1);
  	} else	;
  	D->F=D->shareF[0];
	  
}

void solveCharge(Domain *D,LoadList *LL,double ***Rho,int s,int istart,int iend,double coef,int half,int non_half)
{
  	int i,j,k,ii,jj,kk,maxII,maxJJ,maxKK;
  	int jstart,jend,kstart,kend;
  	double Wx[4],Wy[4],Wz[4],effect1D,effect2D,effect3D;
  	double x,x1,x2,x3,x4,y,y1,y2,y3,y4,z,z1,z2,z3,z4,weight,charge;
  	Particle ***particle;
  	particle=D->particle;
  	ptclList *p;

  	jstart=D->jstart;  jend=D->jend;
  	kstart=D->kstart;  kend=D->kend;

	coef*=LL->density/LL->criticalDensity;

	effect1D=effect2D=effect3D=1.0;
  	 	
   maxII=maxJJ=maxKK=2;
	if(D->dimension==1) {
	 	effect2D=0.0;
	 	effect3D=0.0;
	 	maxJJ=maxKK=1;
	} else if(D->dimension==2) {
	 	effect3D=0.0;
	 	maxKK=1;
	} else ; 	

	for(i=istart; i<iend; i++)
	 	for(j=jstart; j<jend; j++)
	  		for(k=kstart; k<kend; k++)
      	{
        		p=particle[i][j][k].head[s]->pt;
        		while(p)
        		{
          		weight=p->weight;
	      		charge=p->charge;
          		x=half*(p->x+i+p->oldX)*0.5 + non_half*p->x;
					y=half*(p->y+j+p->oldY)*0.5 + non_half*p->y;
					z=half*(p->z+k+p->oldZ)*0.5 + non_half*p->z;
          		Wx[1]=(x-(int)x)*effect1D;	Wx[0]=1-Wx[1];
          		Wy[1]=(y-(int)y)*effect2D;	Wy[0]=1-Wy[1];
          		Wz[1]=(z-(int)z)*effect3D;	Wz[0]=1-Wz[1];
          		for(ii=0; ii<maxII; ii++)
            		for(jj=0; jj<maxJJ; jj++)
              			for(kk=0; kk<maxKK; kk++)
                			Rho[i+ii][j+jj][k+kk]+=Wx[ii]*Wy[jj]*Wz[kk]*coef*weight*charge;
          		p=p->next;
        		}
      	}

}  
