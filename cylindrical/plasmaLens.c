#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_qrng.h>


void plasmaLens(Domain *D,PlasmaLens *PL,int iteration)
{
   int s,n,i,j,numMode,istart,iend,jstart,jend,nSpecies,l;
   int minRSub,minZSub;
   double ne,x,y,z,r,xi,unitXi,Bphi,normalB,m_I,coef,posX,min,max;
   double c[6],xx[6],cosTh,sinTh;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   Particle **particle;
   particle=D->particle;
   LoadList *LL;

   ptclList *p;

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   minRSub=D->minYSub;
   minZSub=D->minXSub;
	nSpecies=D->nSpecies;

   c[0]=0.3155/2.0;
   c[1]=-0.0233967/3.0;
   c[2]=0.0786747/4.0;
   c[3]=-0.878334/5.0;
   c[4]=1.32213/6.0;
   c[5]=-0.751897/7.0;
   unitXi=D->dr*D->lambda/PL->radius;
   normalB=1.0/(eMass*D->omega/(-eCharge));
   m_I=0.106493;
   coef=mu0*PL->current*0.5/M_PI/PL->radius/m_I;

	max=PL->xpoint[PL->xnodes-1];
	min=PL->xpoint[0];
	if( (D->minXDomain>=min && D->minXDomain<max) || (D->maxXDomain>=min && D->maxXDomain<max) )
		istart=istart;
	else 
		istart=iend;

   //position define      
   for(i=istart; i<iend; i++) 
   {
     posX=(double)(i+minZSub-istart);
     for(l=0; l<PL->xnodes-1; l++) {
       if(posX>=PL->xpoint[l] && posX<PL->xpoint[l+1]) 
         ne=((PL->xn[l+1]-PL->xn[l])/(PL->xpoint[l+1]-PL->xpoint[l])*(posX-PL->xpoint[l])+PL->xn[l]);
       else
	      ne=1.0;
     }
     for(j=jstart; j<jend; j++)
         for(s=0; s<nSpecies; s++)  {
           p=particle[i][j].head[s]->pt;
           while(p)  {
			    x=p->x; y=p->y;
             r=sqrt(x*x+y*y); 
	          xi=r*unitXi;
	          xx[0]=xi;
             xx[1]=xx[0]*xi;
             xx[2]=xx[1]*xi;
             xx[3]=xx[2]*xi;
             xx[4]=xx[3]*xi;
             xx[5]=xx[4]*xi;

             Bphi=0.0;
	          for(n=0; n<6; n++) Bphi+=c[n]*xx[n];
	          Bphi*=normalB*coef;
				 cosTh=x/r;
				 sinTh=y/r;
             p->Bx-=sinTh*Bphi*ne;
             p->By+=cosTh*Bphi*ne;

             p=p->next;
	        }
         }
   }
}

