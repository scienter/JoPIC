#include <stdio.h>
#include <stdlib.h>
#include "./include/mesh.h"
#include "./include/constants.h"
#include <math.h>

// PML Left
void settingH_Yee2D_Left(Domain *D,PML *Pml)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    double dt;

    dt=D->dt;    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    k=0;
    for(i=0; i<2; i++)
      for(j=jstart-1; j<=jend; j++)
      {
        Pml->LHx[1-i][j][k]=D->Bx[i+istart][j][k];
      }
    for(i=0; i<1; i++)
      for(j=jstart-1; j<=jend; j++)
      {
        Pml->LHy[i][j][k]=D->By[i+istart][j][k];
        Pml->LHz[i][j][k]=D->Bz[i+istart][j][k];
      }
}

void settingD_Yee2D_Left(Domain *D,PML *Pml)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    double dt;

    dt=D->dt;    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    k=0;
    for(i=0; i<1; i++)
      for(j=jstart-1; j<=jend; j++)
      {
        Pml->LDx[i][j][k]=D->Ex[i+istart][j][k];
      }
    for(i=0; i<2; i++)
      for(j=jstart-1; j<=jend; j++)
      {
        Pml->LDy[1-i][j][k]=D->Ey[i+istart][j][k];
        Pml->LDz[1-i][j][k]=D->Ez[i+istart][j][k];
      }

}


void pmlSolve_DE_Yee2D_Left(Domain *D,PML *Pml)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,pmlCell;
   double dt,dx,dy,coef1,coef2,coef3,freq,alpha,parameter;
   double oldLDx,oldLDy,oldLDz;
   double oldLBx,oldLBy,oldLBz;
   double func2D();

   dt=D->dt;    
   dx=D->dx;    
   dy=D->dy;    
   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   pmlCell=D->pmlCell;
   freq=velocityC/D->lambda;
   parameter=1.0/(pmlCell);

   double sigmaA[pmlCell]; 
   double sigmaB[pmlCell]; 
   for(i=0; i<2; i++) 
     sigmaA[i]=0.0;
   for(i=2; i<pmlCell; i++) 
     sigmaA[i]=func2D(parameter,i,1.0,dx);
   for(i=0; i<1; i++) 
     sigmaB[i]=0.0;
   for(i=1; i<pmlCell; i++) 
     sigmaB[i]=func2D(parameter,i,0.5,dx);
     
   k=0;
   // D, E field solve
   for(i=1; i<pmlCell; i++)
     for(j=jstart; j<jend; j++)
     {
       oldLDx=Pml->LDx[i][j][k];
       oldLDy=Pml->LDy[i][j][k];
       oldLDz=Pml->LDz[i][j][k];

       // D field solve
       coef1=1;
       coef2=dt/dy;
       coef3=1;
       Pml->LDx[i][j][k]=coef1*oldLDx+coef2*(Pml->LHz[i][j][k]-Pml->LHz[i][j-1][k]);
       coef1=1;
       coef2=dt/dx;
       coef3=1;
       Pml->LDy[i][j][k]=coef1*oldLDy-coef2*(Pml->LHz[i][j][k]-Pml->LHz[i-1][j][k]);
       coef1=(2*eps0*freq-sigmaA[i]*dt)/(2*eps0*freq+sigmaA[i]*dt);
       coef2=(2*eps0*freq*dt)/(2*eps0*freq+sigmaA[i]*dt)/dx;
       coef3=(2*eps0*freq*dt)/(2*eps0*freq+sigmaA[i]*dt)/dy;
       Pml->LDz[i][j][k]=coef1*oldLDz+coef2*(Pml->LHy[i][j][k]-Pml->LHy[i-1][j][k])-coef3*(Pml->LHx[i][j][k]-Pml->LHx[i][j-1][k]);

       // E field solve
       coef1=1;
       coef2=(2*eps0*freq+dt*sigmaB[i])/(2*eps0*freq);
       coef3=(2*eps0*freq-dt*sigmaB[i])/(2*eps0*freq);
       Pml->LEx[i][j][k]=coef1*Pml->LEx[i][j][k]+coef2*Pml->LDx[i][j][k]-coef3*oldLDx;
       coef1=(2*eps0*freq-dt*sigmaA[i])/(2*eps0*freq+dt*sigmaA[i]);
       coef2=(2*eps0*freq)/(2*eps0*freq+dt*sigmaA[i]);
       coef3=coef2;
       Pml->LEy[i][j][k]=coef1*Pml->LEy[i][j][k]+coef2*Pml->LDy[i][j][k]-coef3*oldLDy;
       coef1=1;
       coef2=1;
       coef3=1;
       Pml->LEz[i][j][k]=coef1*Pml->LEz[i][j][k]+coef2*Pml->LDz[i][j][k]-coef3*oldLDz;
     }
}

void pmlSolve_BH_Yee2D_Left(Domain *D,PML *Pml)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,pmlCell;
   double dt,dx,dy,coef1,coef2,coef3,freq,alpha,parameter;
   double oldLDx,oldLDy,oldLDz;
   double oldLBx,oldLBy,oldLBz;
   double func2D();

   dt=D->dt;    
   dx=D->dx;    
   dy=D->dy;    
   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   pmlCell=D->pmlCell;
   freq=velocityC/D->lambda;
   parameter=1.0/(pmlCell);

   double sigmaA[pmlCell]; 
   double sigmaB[pmlCell]; 
   for(i=0; i<2; i++) 
     sigmaA[i]=0.0;
   for(i=2; i<pmlCell; i++) 
     sigmaA[i]=func2D(parameter,i,1.0,dx);
   for(i=0; i<1; i++) 
     sigmaB[i]=0.0;
   for(i=1; i<pmlCell; i++) 
     sigmaB[i]=func2D(parameter,i,0.5,dx);
     
   k=0;
   // B, H field solve
   for(i=1; i<pmlCell-1; i++)
     for(j=jstart; j<jend; j++)
     {
       oldLBx=Pml->LBx[i][j][k];
       oldLBy=Pml->LBy[i][j][k];
       oldLBz=Pml->LBz[i][j][k];

       // B field solve
       coef1=1;
       coef2=dt/dy;
       coef3=1;
       Pml->LBx[i][j][k]=coef1*oldLBx-coef2*(Pml->LEz[i][j+1][k]-Pml->LEz[i][j][k]);
       coef1=1;
       coef2=dt/dx;
       coef3=1;
       Pml->LBy[i][j][k]=coef1*oldLBy+coef2*(Pml->LEz[i+1][j][k]-Pml->LEz[i][j][k]);
       coef1=(2*eps0*freq-sigmaB[i]*dt)/(2*eps0*freq+sigmaB[i]*dt);
       coef2=(2*eps0*freq*dt)/(2*eps0*freq+sigmaB[i]*dt)/dx;
       coef3=(2*eps0*freq*dt)/(2*eps0*freq+sigmaB[i]*dt)/dy;
       Pml->LBz[i][j][k]=coef1*oldLBz-coef2*(Pml->LEy[i+1][j][k]-Pml->LEy[i][j][k])+coef3*(Pml->LEx[i][j+1][k]-Pml->LEx[i][j][k]);

       // H field solve
       coef1=1;
       coef2=(2*eps0*freq+dt*sigmaA[i])/(2*eps0*freq);
       coef3=(2*eps0*freq-dt*sigmaA[i])/(2*eps0*freq);
       Pml->LHx[i][j][k]=coef1*Pml->LHx[i][j][k]+coef2*Pml->LBx[i][j][k]-coef3*oldLBx;
       coef1=(2*eps0*freq-dt*sigmaB[i])/(2*eps0*freq+dt*sigmaB[i]);
       coef2=(2*eps0*freq)/(2*eps0*freq+dt*sigmaB[i]);
       coef3=coef2;
       Pml->LHy[i][j][k]=coef1*Pml->LHy[i][j][k]+coef2*Pml->LBy[i][j][k]-coef3*oldLBy;
       coef1=1;
       coef2=1;
       coef3=1;
       Pml->LHz[i][j][k]=coef1*Pml->LHz[i][j][k]+coef2*Pml->LBz[i][j][k]-coef3*oldLBz;
     }
}

void pmlApplyE_Yee2D_Left(Domain *D,PML *Pml)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   k=0;
   for(j=jstart-1; j<=jend; j++)
   {
     D->Ex[istart-1][j][k]=Pml->LEx[1][j][k];
     D->Ey[istart-1][j][k]=Pml->LEy[2][j][k];
     D->Ez[istart-1][j][k]=Pml->LEz[2][j][k];
   }
}

void pmlApplyB_Yee2D_Left(Domain *D,PML *Pml)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   k=0;
   for(j=jstart-1; j<=jend; j++)
   {
     D->Bx[istart-1][j][k]=Pml->LBx[2][j][k];
     D->By[istart-1][j][k]=Pml->LBy[1][j][k];
     D->Bz[istart-1][j][k]=Pml->LBz[1][j][k];
   }
}


// PML Right
void settingH_Yee2D_Right(Domain *D,PML *Pml)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    double dt;

    dt=D->dt;    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    k=0;
   
    for(i=iend-2; i<iend; i++)
      for(j=jstart-1; j<=jend; j++)
      {
        Pml->RHx[i-(iend-2)][j][k]=D->Bx[i][j][k];
      }
    for(i=iend-2; i<iend-1; i++)
      for(j=jstart-1; j<=jend; j++)
      {
        Pml->RHy[i-(iend-2)][j][k]=D->By[i][j][k];
        Pml->RHz[i-(iend-2)][j][k]=D->Bz[i][j][k];
      }
}

void settingD_Yee2D_Right(Domain *D,PML *Pml)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    double dt;

    dt=D->dt;    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    k=0;
    for(i=iend-2; i<iend-1; i++)
      for(j=jstart-1; j<=jend; j++)
      {
        Pml->RDx[i-(iend-2)][j][k]=D->Ex[i][j][k];
      }
    for(i=iend-2; i<iend; i++)
      for(j=jstart-1; j<=jend; j++)
      {
        Pml->RDy[i-(iend-2)][j][k]=D->Ey[i][j][k];
        Pml->RDz[i-(iend-2)][j][k]=D->Ez[i][j][k];
      }
}

void pmlSolve_DE_Yee2D_Right(Domain *D,PML *Pml)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,pmlCell;
   double dt,dx,dy,coef1,coef2,coef3,freq,alpha,parameter;
   double oldRDx,oldRDy,oldRDz;
   double oldRBx,oldRBy,oldRBz;
   double func2D();

   dt=D->dt;    
   dx=D->dx;    
   dy=D->dy;    
   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   pmlCell=D->pmlCell;
   freq=velocityC/D->lambda;
   parameter=1.0/(pmlCell);

   double sigmaA[pmlCell]; 
   double sigmaB[pmlCell]; 
   for(i=0; i<2; i++) 
     sigmaA[i]=0.0;
   for(i=2; i<pmlCell; i++) 
     sigmaA[i]=func2D(parameter,i,1.0,dx);
   for(i=0; i<1; i++) 
     sigmaB[i]=0.0;
   for(i=1; i<pmlCell; i++) 
     sigmaB[i]=func2D(parameter,i,0.5,dx);
     
   k=0;
   // D, E field solve
   for(i=1; i<pmlCell; i++)
     for(j=jstart; j<jend; j++)
     {
       oldRDx=Pml->RDx[i][j][k];
       oldRDy=Pml->RDy[i][j][k];
       oldRDz=Pml->RDz[i][j][k];

       // D field solve
       coef1=1;
       coef2=dt/dy;
       coef3=1;
       Pml->RDx[i][j][k]=coef1*oldRDx+coef2*(Pml->RHz[i][j][k]-Pml->RHz[i][j-1][k]);
       coef1=1;
       coef2=dt/dx;
       coef3=1;
       Pml->RDy[i][j][k]=coef1*oldRDy-coef2*(Pml->RHz[i][j][k]-Pml->RHz[i-1][j][k]);
       coef1=(2*eps0*freq-sigmaA[i]*dt)/(2*eps0*freq+sigmaA[i]*dt);
       coef2=(2*eps0*freq*dt)/(2*eps0*freq+sigmaA[i]*dt)/dx;
       coef3=(2*eps0*freq*dt)/(2*eps0*freq+sigmaA[i]*dt)/dy;
       Pml->RDz[i][j][k]=coef1*oldRDz+coef2*(Pml->RHy[i][j][k]-Pml->RHy[i-1][j][k])-coef3*(Pml->RHx[i][j][k]-Pml->RHx[i][j-1][k]);

       // E field solve
       coef1=1;
       coef2=(2*eps0*freq+dt*sigmaB[i])/(2*eps0*freq);
       coef3=(2*eps0*freq-dt*sigmaB[i])/(2*eps0*freq);
       Pml->REx[i][j][k]=coef1*Pml->REx[i][j][k]+coef2*Pml->RDx[i][j][k]-coef3*oldRDx;
       coef1=(2*eps0*freq-dt*sigmaA[i])/(2*eps0*freq+dt*sigmaA[i]);
       coef2=(2*eps0*freq)/(2*eps0*freq+dt*sigmaA[i]);
       coef3=coef2;
       Pml->REy[i][j][k]=coef1*Pml->REy[i][j][k]+coef2*Pml->RDy[i][j][k]-coef3*oldRDy;
       coef1=1;
       coef2=1;
       coef3=1;
       Pml->REz[i][j][k]=coef1*Pml->REz[i][j][k]+coef2*Pml->RDz[i][j][k]-coef3*oldRDz;
     }
}

void pmlSolve_BH_Yee2D_Right(Domain *D,PML *Pml)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,pmlCell;
   double dt,dx,dy,coef1,coef2,coef3,freq,alpha,parameter;
   double oldRDx,oldRDy,oldRDz;
   double oldRBx,oldRBy,oldRBz;
   double func2D();

   dt=D->dt;    
   dx=D->dx;    
   dy=D->dy;    
   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   pmlCell=D->pmlCell;
   freq=velocityC/D->lambda;
   parameter=1.0/(pmlCell);

   double sigmaA[pmlCell]; 
   double sigmaB[pmlCell]; 
   for(i=0; i<2; i++) 
     sigmaA[i]=0.0;
   for(i=2; i<pmlCell; i++) 
     sigmaA[i]=func2D(parameter,i,1.0,dx);
   for(i=0; i<1; i++) 
     sigmaB[i]=0.0;
   for(i=1; i<pmlCell; i++) 
     sigmaB[i]=func2D(parameter,i,0.5,dx);
     
   k=0;

   // B, H field solve
   for(i=1; i<pmlCell-1; i++)
     for(j=jstart; j<jend; j++)
     {
       oldRBx=Pml->RBx[i][j][k];
       oldRBy=Pml->RBy[i][j][k];
       oldRBz=Pml->RBz[i][j][k];

       // B field solve
       coef1=1;
       coef2=dt/dy;
       coef3=1;
       Pml->RBx[i][j][k]=coef1*oldRBx-coef2*(Pml->REz[i][j+1][k]-Pml->REz[i][j][k]);
       coef1=1;
       coef2=dt/dx;
       coef3=1;
       Pml->RBy[i][j][k]=coef1*oldRBy+coef2*(Pml->REz[i+1][j][k]-Pml->REz[i][j][k]);
       coef1=(2*eps0*freq-sigmaB[i]*dt)/(2*eps0*freq+sigmaB[i]*dt);
       coef2=(2*eps0*freq*dt)/(2*eps0*freq+sigmaB[i]*dt)/dx;
       coef3=(2*eps0*freq*dt)/(2*eps0*freq+sigmaB[i]*dt)/dy;
       Pml->RBz[i][j][k]=coef1*oldRBz-coef2*(Pml->REy[i+1][j][k]-Pml->REy[i][j][k])+coef3*(Pml->REx[i][j+1][k]-Pml->REx[i][j][k]);

       // H field solve
       coef1=1;
       coef2=(2*eps0*freq+dt*sigmaA[i])/(2*eps0*freq);
       coef3=(2*eps0*freq-dt*sigmaA[i])/(2*eps0*freq);
       Pml->RHx[i][j][k]=coef1*Pml->RHx[i][j][k]+coef2*Pml->RBx[i][j][k]-coef3*oldRBx;
       coef1=(2*eps0*freq-dt*sigmaB[i])/(2*eps0*freq+dt*sigmaB[i]);
       coef2=(2*eps0*freq)/(2*eps0*freq+dt*sigmaB[i]);
       coef3=coef2;
       Pml->RHy[i][j][k]=coef1*Pml->RHy[i][j][k]+coef2*Pml->RBy[i][j][k]-coef3*oldRBy;
       coef1=1;
       coef2=1;
       coef3=1;
       Pml->RHz[i][j][k]=coef1*Pml->RHz[i][j][k]+coef2*Pml->RBz[i][j][k]-coef3*oldRBz;
     }
}

void pmlApplyE_Yee2D_Right(Domain *D,PML *Pml)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   k=0;
   for(j=jstart-1; j<=jend; j++)
   {
     D->Ex[iend-1][j][k]=Pml->REx[1][j][k];
     D->Ey[iend][j][k]=Pml->REy[2][j][k];
     D->Ez[iend][j][k]=Pml->REz[2][j][k];
   }
}

void pmlApplyB_Yee2D_Right(Domain *D,PML *Pml)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   k=0;
   for(j=jstart-1; j<=jend; j++)
   {
     D->Bx[iend][j][k]=Pml->RBx[2][j][k];
     D->By[iend-1][j][k]=Pml->RBy[1][j][k];
     D->Bz[iend-1][j][k]=Pml->RBz[1][j][k];
   }
}

double func2D(double parameter, int i,double offset, double dx)
{
  double tmp,result;

  tmp=(i-offset)*parameter;
  result=pow(tmp,2)*300000;
     
  return result;
}
