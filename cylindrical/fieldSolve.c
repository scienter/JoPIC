#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "mesh.h"
#include "constants.h"
#include "math.h"

void absorb_U(Domain *D);
void Bsolve_Yee(Domain *D,int iteration);
void Esolve_Yee(Domain *D,int iteration);
void Bsolve_NoCherenkov(Domain *D,int iteration);
void EzBz_solve_Split(Domain *D,int iteration);
void PS_solve_Split(Domain *D,int iteration);
void filter_center(Domain *D,int numField,int iteration);
void EzBz_Now_Split(Domain *D,int iteration);
void Yee_Now(Domain *D,int iteration);

void fieldSolve2(Domain D,double t,int iteration)
{
  	int rankX,rankY;
  	float limit;
  	LaserList *L;
  	int myrank, nTasks,rank,rankM,rankN;

  	MPI_Status status;
  	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  	switch(D.fieldType) {

  	case Yee :
      //load laser
    	if(D.boostOn==OFF)   {
    	   L=D.laserList;
    	   while(L->next)  {
    	      limit=((L->rU+L->rD)/D.dtRatio*2.1+L->retard)*1.0/D.divisionLambda;
    	      if(iteration<=limit && L->loadMethod==Boundary) loadLaser(&D,L,t); else ;
    	      L=L->next;
    	   }
    	} else ;

    	Bsolve_Yee(&D,iteration);
   	D.shareF[0]=D.BrR;
   	D.shareF[1]=D.BrI;
   	D.shareF[2]=D.BpR;
   	D.shareF[3]=D.BpI;
   	D.shareF[4]=D.BzR;
   	D.shareF[5]=D.BzI;
   	MPI_TransferFNew_Xplus(&D,6,D.nySub+5,3);
   	MPI_TransferFNew_Xminus(&D,6,D.nySub+5,3);
	 	if(D.Period==ON) { MPI_TransferFNew_Period_X(&D,6,D.nySub+5,3); } else ;  
   	D.BrR=D.shareF[0];
   	D.BrI=D.shareF[1];
   	D.BpR=D.shareF[2];
   	D.BpI=D.shareF[3];
   	D.BzR=D.shareF[4];
   	D.BzI=D.shareF[5];

   	D.shareF[0]=D.BrR;
   	D.shareF[1]=D.BrI;
   	D.shareF[2]=D.BpR;
   	D.shareF[3]=D.BpI;
   	D.shareF[4]=D.BzR;
   	D.shareF[5]=D.BzI;
      filter_center(&D,6,iteration);
      MPI_TransferFNew_Xplus(&D,6,D.nySub+5,3);
      MPI_TransferFNew_Xminus(&D,6,D.nySub+5,3);
      if(D.Period==ON) { MPI_TransferFNew_Period_X(&D,6,D.nySub+5,3); } else ;
   	D.BrR=D.shareF[0];
   	D.BrI=D.shareF[1];
   	D.BpR=D.shareF[2];
   	D.BpI=D.shareF[3];
   	D.BzR=D.shareF[4];
   	D.BzI=D.shareF[5];
    	//MPI_Transfer12F_Xminus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
    	//MPI_Transfer12F_Xplus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
    	//if(D.Period==ON)
    	//  MPI_Transfer12F_Period_X(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
      
      Yee_Now(&D,iteration);
    	break ;

	case NoCherenkov :
    	//load laser
		  if(D.boostOn==OFF)   {
				L=D.laserList;
				while(L->next)  {
					 limit=((L->rU+L->rD)/D.dtRatio*2.1+L->retard)*1.0/D.divisionLambda;
					 if(iteration<=limit && L->loadMethod==Boundary) loadLaser(&D,L,t); else ;
					 L=L->next;
				}
		  } else ;

		Bsolve_NoCherenkov(&D,iteration);
		MPI_Transfer12F_Xminus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
		MPI_Transfer12F_Xplus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
		if(D.Period==ON)
		MPI_Transfer12F_Period_X(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
		else ;

	   break ;

	case Split :
		EzBz_solve_Split(&D,iteration);
   	D.shareF[0]=D.EzR;
   	D.shareF[1]=D.EzI;
   	D.shareF[2]=D.BzR;
   	D.shareF[3]=D.BzI;
   	MPI_TransferFNew_Xplus(&D,4,D.nySub+5,3);
   	MPI_TransferFNew_Xminus(&D,4,D.nySub+5,3);
	 	if(D.Period==ON) { MPI_TransferFNew_Period_X(&D,4,D.nySub+5,3); } else ;  
   	D.EzR=D.shareF[0];
   	D.EzI=D.shareF[1];
   	D.BzR=D.shareF[2];
   	D.BzI=D.shareF[3];

   	D.shareF[0]=D.EzR;
   	D.shareF[1]=D.EzI;
   	D.shareF[2]=D.BzR;
   	D.shareF[3]=D.BzI;
      filter_center(&D,4,iteration);
      MPI_TransferFNew_Xplus(&D,4,D.nySub+5,3);
      MPI_TransferFNew_Xminus(&D,4,D.nySub+5,3);
      if(D.Period==ON) { MPI_TransferFNew_Period_X(&D,4,D.nySub+5,3); } else ;
      D.EzR=D.shareF[0];
      D.EzI=D.shareF[1];
      D.BzR=D.shareF[2];
      D.BzI=D.shareF[3];      

      EzBz_Now_Split(&D,iteration);
      
		break ;
	}
}

void fieldSolve1(Domain D,double t,int iteration)
{
  	int rankX,rankY,numField;
  	float limit;
  	LaserList *L;
  	int myrank, nTasks,rank,rankM,rankN;

  	MPI_Status status;
  	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   //load laser
   if(D.boostOn==OFF)   {
     L=D.laserList;
     while(L->next)  {
       limit=((L->rU+L->rD)/D.dtRatio*2.1+L->retard)*1.0/D.divisionLambda;
       if(iteration<=limit && L->loadMethod==Boundary) loadLaser(&D,L,t); else ;
       L=L->next;
     }
   } else ;

  	switch(D.fieldType) {
  	case Yee :
  	  	Esolve_Yee(&D,iteration);
   	D.shareF[0]=D.ErR;
   	D.shareF[1]=D.ErI;
   	D.shareF[2]=D.EpR;
   	D.shareF[3]=D.EpI;
   	D.shareF[4]=D.EzR;
   	D.shareF[5]=D.EzI;
   	MPI_TransferFNew_Xplus(&D,6,D.nySub+5,3);
   	MPI_TransferFNew_Xminus(&D,6,D.nySub+5,3);
	 	if(D.Period==ON) { MPI_TransferFNew_Period_X(&D,6,D.nySub+5,3); } else ;  
   	D.ErR=D.shareF[0];
   	D.ErI=D.shareF[1];
   	D.EpR=D.shareF[2];
   	D.EpI=D.shareF[3];
   	D.EzR=D.shareF[4];
   	D.EzI=D.shareF[5];

   	D.shareF[0]=D.ErR;
   	D.shareF[1]=D.ErI;
   	D.shareF[2]=D.EpR;
   	D.shareF[3]=D.EpI;
   	D.shareF[4]=D.EzR;
   	D.shareF[5]=D.EzI;
      filter_center(&D,6,iteration);
      MPI_TransferFNew_Xplus(&D,6,D.nySub+5,3);
      MPI_TransferFNew_Xminus(&D,6,D.nySub+5,3);
      if(D.Period==ON) { MPI_TransferFNew_Period_X(&D,6,D.nySub+5,3); } else ;
   	D.ErR=D.shareF[0];
   	D.ErI=D.shareF[1];
   	D.EpR=D.shareF[2];
   	D.EpI=D.shareF[3];
   	D.EzR=D.shareF[4];
   	D.EzI=D.shareF[5];
  	  	break ;

  	case NoCherenkov :
    Esolve_Yee(&D,iteration);
    MPI_Transfer6F_Xminus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
    MPI_Transfer6F_Xplus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
    if(D.Period==ON)
      MPI_Transfer6F_Period_X(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
	 else ;
    break ;

  	case Split :
   	PS_solve_Split(&D,iteration);
   	D.shareF[0]=D.PrR;
   	D.shareF[1]=D.PrI;
   	D.shareF[2]=D.SrR;
   	D.shareF[3]=D.SrI;
   	D.shareF[4]=D.PlR;
   	D.shareF[5]=D.PlI;
   	D.shareF[6]=D.SlR;
   	D.shareF[7]=D.SlI;
   	MPI_TransferFNew_Xplus(&D,8,D.nySub+5,3);
   	MPI_TransferFNew_Xminus(&D,8,D.nySub+5,3);
		 if(D.Period==ON) { MPI_TransferFNew_Period_X(&D,8,D.nySub+5,3); } else ;  
   	D.PrR=D.shareF[0];
   	D.PrI=D.shareF[1];
   	D.SrR=D.shareF[2];
   	D.SrI=D.shareF[3];
   	D.PlR=D.shareF[4];
   	D.PlI=D.shareF[5];
   	D.SlR=D.shareF[6];
   	D.SlI=D.shareF[7]; 

   	D.shareF[0]=D.PrR;
   	D.shareF[1]=D.PrI;
   	D.shareF[2]=D.SrR;
   	D.shareF[3]=D.SrI;
   	D.shareF[4]=D.PlR;
   	D.shareF[5]=D.PlI;
   	D.shareF[6]=D.SlR;
   	D.shareF[7]=D.SlI;
         filter_center(&D,8,iteration);
         MPI_TransferFNew_Xplus(&D,8,D.nySub+5,3);
         MPI_TransferFNew_Xminus(&D,8,D.nySub+5,3);
         if(D.Period==ON) { MPI_TransferFNew_Period_X(&D,8,D.nySub+5,3); } else ;
   	D.PrR=D.shareF[0];
   	D.PrI=D.shareF[1];
   	D.SrR=D.shareF[2];
   	D.SrI=D.shareF[3];
   	D.PlR=D.shareF[4];
   	D.PlI=D.shareF[5];
   	D.SlR=D.shareF[6];
   	D.SlI=D.shareF[7]; 

   	break ;
  	}
}

void Bsolve_NoCherenkov(Domain *D,int iteration)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double x,minZSub,dtBydr,dtBydz,r,dr,dz,dt,BzR,BrR,BpR,BzI,BrI,BpI;
  double oldBzR,oldBrR,oldBpR,oldBzI,oldBrI,oldBpI;
  double deltaZ,betaPZ,betaRZ,alphaZ,alphaP,alphaR;
  double upr,upd,lftr,lftd,upL,leftL,LdU,LdL,rr,rd,tmp,tmpr,tmpd;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;
  minZSub=D->minXSub;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr; dz=D->dz; dt=D->dt;

  leftL=(double)(D->minXDomain+D->pmlCellLeft);
  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  LdL=D->pmlCellLeft;
  rr=D->pmlr;    rd=D->pmld;
  lftr=lftd=upr=upd=1.0;

  tmp=sin(0.5*pi*dtBydz);
  deltaZ=0.25*(1.0-tmp*tmp/dtBydz/dtBydz);
  betaPZ=0.25;
  betaRZ=0.25;
  alphaZ=1.0-3.0*deltaZ;
  alphaP=1.0-2.0*betaPZ;
  alphaR=1.0-2.0*betaRZ;

  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)
    {
      x=(i-istart)+minZSub;
      for(j=jstart+1; j<jend; j++)  
      {
        r=(double)(j-jstart);
        oldBzR=D->BzR[m][i][j];
        oldBrR=D->BrR[m][i][j];
        oldBpR=D->BpR[m][i][j];
        oldBzI=D->BzI[m][i][j];
        oldBrI=D->BrI[m][i][j];
        oldBpI=D->BpI[m][i][j];

        tmpr=upr*lftr;
        tmpd=upd*lftd;
        tmp=tmpr*dtBydr/r*(-alphaR*((r+0.5)*D->EpR[m][i][j+1]-(r-0.5)*D->EpR[m][i][j])-betaRZ*((r+0.5)*D->EpR[m][i+1][j+1]-(r-0.5)*D->EpR[m][i+1][j]+(r+0.5)*D->EpR[m][i-1][j+1]-(r-0.5)*D->EpR[m][i-1][j])-m*(alphaP*D->ErI[m][i][j]+betaPZ*(D->ErI[m][i+1][j]+D->ErI[m][i-1][j])));
        D->BzR[m][i][j]=tmpd*(oldBzR+tmp);
        tmp=tmpr*dtBydr/r*(-alphaR*((r+0.5)*D->EpI[m][i][j+1]-(r-0.5)*D->EpI[m][i][j])-betaRZ*((r+0.5)*D->EpI[m][i+1][j+1]-(r-0.5)*D->EpI[m][i+1][j]+(r+0.5)*D->EpI[m][i-1][j+1]-(r-0.5)*D->EpI[m][i-1][j])+m*(alphaP*D->ErR[m][i][j]+betaPZ*(D->ErR[m][i+1][j]+D->ErR[m][i-1][j])));
        D->BzI[m][i][j]=tmpd*(oldBzI+tmp);

        tmp=tmpr*(dtBydz*(alphaZ*(D->EpR[m][i+1][j]-D->EpR[m][i][j])+deltaZ*(D->EpR[m][i+2][j]-D->EpR[m][i-1][j]))+m*dtBydr/(r-0.5)*(alphaP*D->EzI[m][i][j]+betaPZ*(D->EzI[m][i+1][j]+D->EzI[m][i-1][j])));
        D->BrR[m][i][j]=tmpd*(oldBrR+tmp);
        tmp=tmpr*(dtBydz*(alphaZ*(D->EpI[m][i+1][j]-D->EpI[m][i][j])+deltaZ*(D->EpI[m][i+2][j]-D->EpI[m][i-1][j]))-m*dtBydr/(r-0.5)*(alphaP*D->EzR[m][i][j]+betaPZ*(D->EzR[m][i+1][j]+D->EzR[m][i-1][j])));
        D->BrI[m][i][j]=tmpd*(oldBrI+tmp);

        tmp=tmpr*(dtBydr*(alphaR*(D->EzR[m][i][j+1]-D->EzR[m][i][j])+betaRZ*(D->EzR[m][i+1][j+1]-D->EzR[m][i+1][j]+D->EzR[m][i-1][j+1]-D->EzR[m][i-1][j]))-dtBydz*(alphaZ*(D->ErR[m][i+1][j]-D->ErR[m][i][j])+deltaZ*(D->ErR[m][i+2][j]-D->ErR[m][i-1][j])));
        D->BpR[m][i][j]=tmpd*(oldBpR+tmp);
        tmp=tmpr*(dtBydr*(alphaR*(D->EzI[m][i][j+1]-D->EzI[m][i][j])+betaRZ*(D->EzI[m][i+1][j+1]-D->EzI[m][i+1][j]+D->EzI[m][i-1][j+1]-D->EzI[m][i-1][j]))-dtBydz*(alphaZ*(D->ErI[m][i+1][j]-D->ErI[m][i][j])+deltaZ*(D->ErI[m][i+2][j]-D->ErI[m][i-1][j])));
        D->BpI[m][i][j]=tmpd*(oldBpI+tmp);

        D->BzNowR[m][i][j]=0.5*(D->BzR[m][i][j]+oldBzR);
        D->BrNowR[m][i][j]=0.5*(D->BrR[m][i][j]+oldBrR);
        D->BpNowR[m][i][j]=0.5*(D->BpR[m][i][j]+oldBpR);
        D->BzNowI[m][i][j]=0.5*(D->BzI[m][i][j]+oldBzI);
        D->BrNowI[m][i][j]=0.5*(D->BrI[m][i][j]+oldBrI);
        D->BpNowI[m][i][j]=0.5*(D->BpI[m][i][j]+oldBpI);
      }
    }

  lftr=lftd=upr=upd=1.0;
  j=jstart;
    for(i=istart; i<iend; i++) {
      x=(i-istart)+minZSub;

      oldBzR=D->BzR[0][i][j];
      oldBpR=D->BpR[1][i][j];
      oldBpI=D->BpI[1][i][j];

      tmpr=upr*lftr;
      tmpd=upd*lftd;

      D->BzR[0][i][j]+=-4.0*dtBydr*D->EpR[0][i][j+1];

      tmp=tmpr*(dt*(4.0*D->EzR[1][i][j+1]-D->EzR[1][i][j+2])-dtBydz*(alphaZ*(D->ErR[1][i+1][j]-D->ErR[1][i][j])+deltaZ*(D->ErR[1][i+2][j]-D->ErR[1][i-1][j])));
      D->BpR[1][i][j]=tmpd*(oldBpR+tmp);
      tmp=tmpr*(dt*(4.0*D->EzI[1][i][j+1]-D->EzI[1][i][j+2])-dtBydz*(alphaZ*(D->ErI[1][i+1][j]-D->ErI[1][i][j])+deltaZ*(D->ErI[1][i+2][j]-D->ErI[1][i-1][j])));
      D->BpI[1][i][j]=tmpd*(oldBpI+tmp);
      D->BzNowR[0][i][j]=0.5*(D->BzR[0][i][j]+oldBzR);
      D->BpNowR[1][i][j]=0.5*(D->BpR[1][i][j]+oldBpR);
      D->BpNowI[1][i][j]=0.5*(D->BpI[1][i][j]+oldBpI);
    }
}

void Bsolve_Yee(Domain *D,int iteration)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,BzR,BrR,BpR,BzI,BrI,BpI;
  double oldBzR,oldBrR,oldBpR,oldBzI,oldBrI,oldBpI;
  double upr,upd,upL,LdU,rr,rd,tmp,tmpr,tmpd,coef1,coef2,coef3,coef4;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  rr=D->pmlr;    rd=D->pmld;
  upr=upd=1.0;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr;
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)
    {
      for(j=jstart+1; j<jend; j++)  
      {
        r=(double)(j-jstart);
        D->BzNowR[m][i][j]=oldBzR=D->BzR[m][i][j];
        D->BrNowR[m][i][j]=oldBrR=D->BrR[m][i][j];
        D->BpNowR[m][i][j]=oldBpR=D->BpR[m][i][j];
        D->BzNowI[m][i][j]=oldBzI=D->BzI[m][i][j];
        D->BrNowI[m][i][j]=oldBrI=D->BrI[m][i][j];
        D->BpNowI[m][i][j]=oldBpI=D->BpI[m][i][j];

        tmpr=upr;    tmpd=upd;

        tmp=tmpr*(-dtBydr*m/(r+0.5)*D->ErI[m][i][j]
             -0.5*dtBydr/(r+0.5)*(D->EpR[m][i][j+1]+D->EpR[m][i][j])
             -dtBydr*(D->EpR[m][i][j+1]-D->EpR[m][i][j]));
        D->BzR[m][i][j]=tmpd*(oldBzR+tmp);

        tmp=tmpr*(dtBydr*m/(r+0.5)*D->ErR[m][i][j]
             -0.5*dtBydr/(r+0.5)*(D->EpI[m][i][j+1]+D->EpI[m][i][j])
             -dtBydr*(D->EpI[m][i][j+1]-D->EpI[m][i][j]));
        D->BzI[m][i][j]=tmpd*(oldBzI+tmp);

        tmp=tmpr*(dtBydz*(D->EpR[m][i+1][j]-D->EpR[m][i][j])
              +dtBydr*m/r*D->EzI[m][i][j]);
        D->BrR[m][i][j]=tmpd*(oldBrR+tmp);
        tmp=tmpr*(dtBydz*(D->EpI[m][i+1][j]-D->EpI[m][i][j])
              -dtBydr*m/r*D->EzR[m][i][j]);
        D->BrI[m][i][j]=tmpd*(oldBrI+tmp);

        tmp=tmpr*(dtBydr*(D->EzR[m][i][j+1]-D->EzR[m][i][j])
             -dtBydz*(D->ErR[m][i+1][j]-D->ErR[m][i][j]));
        D->BpR[m][i][j]=tmpd*(oldBpR+tmp);
        tmp=tmpr*(dtBydr*(D->EzI[m][i][j+1]-D->EzI[m][i][j])
             -dtBydz*(D->ErI[m][i+1][j]-D->ErI[m][i][j]));
        D->BpI[m][i][j]=tmpd*(oldBpI+tmp);

      }
    }

  j=jstart; m=1;
    for(i=istart; i<iend; i++) {
      D->BrNowR[m][i][j]=oldBrR=D->BrR[m][i][j];
      D->BrNowI[m][i][j]=oldBrI=D->BrI[m][i][j];

      D->BrR[m][i][j]+=dtBydz*(D->EpR[m][i+1][j]-D->EpR[m][i][j])
        +dtBydr*D->EzI[m][i][j+1];
      D->BrI[m][i][j]+=dtBydz*(D->EpI[m][i+1][j]-D->EpI[m][i][j])
        -dtBydr*D->EzR[m][i][j+1];
    }
  r=0.0;
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)
    {
        D->BzNowR[m][i][j]=oldBzR=D->BzR[m][i][j];
        D->BpNowR[m][i][j]=oldBpR=D->BpR[m][i][j];
        D->BzNowI[m][i][j]=oldBzI=D->BzI[m][i][j];
        D->BpNowI[m][i][j]=oldBpI=D->BpI[m][i][j];

        D->BzR[m][i][j]+=-dtBydr*m/(r+0.5)*D->ErI[m][i][j]
          -0.5*dtBydr/(r+0.5)*(D->EpR[m][i][j+1]+D->EpR[m][i][j])
          -dtBydr*(D->EpR[m][i][j+1]-D->EpR[m][i][j]);
        D->BzI[m][i][j]+=dtBydr*m/(r+0.5)*D->ErR[m][i][j]
          -0.5*dtBydr/(r+0.5)*(D->EpI[m][i][j+1]+D->EpI[m][i][j])
          -dtBydr*(D->EpI[m][i][j+1]-D->EpI[m][i][j]);

        D->BpR[m][i][j]+=dtBydr*(D->EzR[m][i][j+1]-D->EzR[m][i][j])
          -dtBydz*(D->ErR[m][i+1][j]-D->ErR[m][i][j]);
        D->BpI[m][i][j]+=dtBydr*(D->EzI[m][i][j+1]-D->EzI[m][i][j])
          -dtBydz*(D->ErI[m][i+1][j]-D->ErI[m][i][j]);
    }

}

void Esolve_Yee(Domain *D,int iteration)
{
  int i,j,m,numMode,istart,iend,jstart,jend,a;  
  double dtBydr,dtBydz,r,dt,dr,dz,dF;
  double EzR,ErR,EpR,EzI,ErI,EpI,x,minZSub;
  double upr,upd,lftr,lftd,upL,leftL,LdU,LdL,rr,rd,tmp,tmpr,tmpd,alpha[2];
  double oldEzR,oldEzI,oldErR,oldErI,oldEpR,oldEpI;
  double beforeEpR,beforeEpI,beforeEzR,beforeEzI,beforeErR,beforeErI;
  double nowEpR,nowEpI,nowEzR,nowEzI,nowErR,nowErI;
  double beforeBpR,beforeBpI,beforeBzR,beforeBzI,beforeBrR,beforeBrI;
  double nowBpR,nowBpI,nowBzR,nowBzI,nowBrR,nowBrI;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;
  minZSub=D->minXSub;

  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	 dF=D->dF;
  leftL=(double)(D->minXDomain+D->pmlCellLeft);
  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  LdL=D->pmlCellLeft;
  rr=D->pmlr;    rd=D->pmld;
  lftr=lftd=upr=upd=1.0;

  dtBydr=D->dt/D->dr; dr=D->dr; dz=D->dz;
  dtBydz=D->dt/D->dz; dt=D->dt;
  m=0;
    for(i=istart; i<iend; i++) 
    {
      for(j=jstart+1; j<jend; j++) 
      {
        r=(double)(j-jstart);

        oldEzR=D->EzR[m][i][j];
        oldEzI=D->EzI[m][i][j];
        oldErR=D->ErR[m][i][j];
        oldErI=D->ErI[m][i][j];
        oldEpR=D->EpR[m][i][j];
        oldEpI=D->EpI[m][i][j];

        tmpr=upr*lftr;     tmpd=upd*lftd;

        tmp=tmpr*(
              0.5*dtBydr/r*(D->BpR[m][i][j]+D->BpR[m][i][j-1])
              +dtBydr*(D->BpR[m][i][j]-D->BpR[m][i][j-1])
              -2.0*pi*dt*D->JzR[m][i][j]
              +dF*dtBydz*(D->FR[m][i+1][j]-D->FR[m][i][j]) 
              );
        D->EzR[m][i][j]=tmpd*(oldEzR+tmp);

        tmp=tmpr*(
              -2.0*pi*dt*D->JrR[m][i][j]
              -dtBydz*(D->BpR[m][i][j]-D->BpR[m][i-1][j])
              +dF*dtBydr*(D->FR[m][i][j+1]-D->FR[m][i][j])
              );
        D->ErR[m][i][j]=tmpd*(oldErR+tmp);
        
        tmp=tmpr*(
              -2.0*pi*dt*D->JpR[m][i][j]
              -dtBydr*(D->BzR[m][i][j]-D->BzR[m][i][j-1])
              +dtBydz*(D->BrR[m][i][j]-D->BrR[m][i-1][j])
              );
        D->EpR[m][i][j]=tmpd*(oldEpR+tmp);
      }
    }
  for(m=1; m<numMode; m++) 
    for(i=istart; i<iend; i++) 
    {
      for(j=jstart+1; j<jend; j++) 
      {
        r=(double)(j-jstart);

        oldEzR=D->EzR[m][i][j];
        oldEzI=D->EzI[m][i][j];
        oldErR=D->ErR[m][i][j];
        oldErI=D->ErI[m][i][j];
        oldEpR=D->EpR[m][i][j];
        oldEpI=D->EpI[m][i][j];

        tmpr=upr*lftr;     tmpd=upd*lftd;

        tmp=tmpr*(
              dtBydr*m/r*D->BrI[m][i][j]
              +0.5*dtBydr/r*(D->BpR[m][i][j]+D->BpR[m][i][j-1])
              +dtBydr*(D->BpR[m][i][j]-D->BpR[m][i][j-1])
              -2.0*pi*dt*D->JzR[m][i][j]
              +dF*dtBydz*(D->FR[m][i+1][j]-D->FR[m][i][j])
              );
        D->EzR[m][i][j]=tmpd*(oldEzR+tmp);
        tmp=tmpr*(
              -dtBydr*m/r*D->BrR[m][i][j]
              +0.5*dtBydr/r*(D->BpI[m][i][j]+D->BpI[m][i][j-1])
              +dtBydr*(D->BpI[m][i][j]-D->BpI[m][i][j-1])
              -2.0*pi*dt*D->JzI[m][i][j]
              +dF*dtBydz*(D->FI[m][i+1][j]-D->FI[m][i][j])
              );
        D->EzI[m][i][j]=tmpd*(oldEzI+tmp);
       
        tmp=tmpr*(
              -2.0*pi*dt*D->JrR[m][i][j]
              -dtBydz*(D->BpR[m][i][j]-D->BpR[m][i-1][j])
              -dtBydr*m/(r+0.5)*D->BzI[m][i][j]
              +dF*dtBydr*(D->FR[m][i][j+1]-D->FR[m][i][j])
              );
        D->ErR[m][i][j]=tmpd*(oldErR+tmp);
        tmp=tmpr*(
              -2.0*pi*dt*D->JrI[m][i][j]
              -dtBydz*(D->BpI[m][i][j]-D->BpI[m][i-1][j])
              +dtBydr*m/(r+0.5)*D->BzR[m][i][j]
              +dF*dtBydr*(D->FI[m][i][j+1]-D->FI[m][i][j])
              );
        D->ErI[m][i][j]=tmpd*(oldErI+tmp);
        
        tmp=tmpr*(
              -2.0*pi*dt*D->JpR[m][i][j]
              -dtBydr*(D->BzR[m][i][j]-D->BzR[m][i][j-1])
              +dtBydz*(D->BrR[m][i][j]-D->BrR[m][i-1][j])
              -dF*dtBydr/r*m*D->FI[m][i][j]
              );
        D->EpR[m][i][j]=tmpd*(oldEpR+tmp);
        tmp=tmpr*(
              -2.0*pi*dt*D->JpI[m][i][j]
              -dtBydr*(D->BzI[m][i][j]-D->BzI[m][i][j-1])
              +dtBydz*(D->BrI[m][i][j]-D->BrI[m][i-1][j])
              +dF*dtBydr/r*m*D->FR[m][i][j]
              );
        D->EpI[m][i][j]=tmpd*(oldEpI+tmp);
      }
    }

  j=jstart;
    for(i=istart; i<iend; i++) {
      m=1;
      D->EpR[m][i][j]+=-2.0*pi*dt*D->JpR[m][i][j]
        -2.0*dtBydr*(D->BzR[m][i][j])
        +dtBydz*(D->BrR[m][i][j]-D->BrR[m][i-1][j])
        -dF*dtBydr*m*D->FI[m][i][j+1];
      D->EpI[m][i][j]+=-2.0*pi*dt*D->JpI[m][i][j]
        -2.0*dtBydr*(D->BzI[m][i][j])
        +dtBydz*(D->BrI[m][i][j]-D->BrI[m][i-1][j])
        +dF*dtBydr*m*D->FR[m][i][j+1];
      m=0;
      D->EzR[m][i][j]+=4.0*dtBydr*D->BpR[m][i][j]
        -2.0*pi*dt*D->JzR[m][i][j]
        +dF*dtBydz*(D->FR[m][i+1][j]-D->FR[m][i][j]);
    }
  for(m=0; m<numMode; m++) 
    for(i=istart; i<iend; i++) 
    {
        r=(double)(j-jstart);
        oldErR=D->ErR[m][i][j];
        oldErI=D->ErI[m][i][j];

        D->ErR[m][i][j]+=-2.0*pi*dt*D->JrR[m][i][j]
          -dtBydz*(D->BpR[m][i][j]-D->BpR[m][i-1][j])
          -dtBydr*m/(r+0.5)*D->BzI[m][i][j]
          +dF*dtBydr*(D->FR[m][i][j+1]-D->FR[m][i][j]);
        D->ErI[m][i][j]+=-2.0*pi*dt*D->JrI[m][i][j]
          -dtBydz*(D->BpI[m][i][j]-D->BpI[m][i-1][j])
          +dtBydr*m/(r+0.5)*D->BzR[m][i][j] 
          +dF*dtBydr*(D->FI[m][i][j+1]-D->FI[m][i][j]);
    }

}

void EzBz_solve_Split(Domain *D,int iteration)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,dz,dt,dF;
  double upPrR,upPlR,upSrR,upSlR;
  double dnPrR,dnPlR,dnSrR,dnSlR,PrR,PlR,SrR,SlR;
  double upPrI,upPlI,upSrI,upSlI;
  double dnPrI,dnPlI,dnSrI,dnSlI,PrI,PlI,SrI,SlI;
  double rr,rd,upL,LdU,tmpr,tmpd,tmp,upr,upd;
  double BzROld,BzIOld,EzROld,EzIOld;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr; dz=D->dz; dt=D->dt; dF=D->dF;

  for(j=jstart; j<jend; j++) {
    D->upr[j]=1.0;
    D->upd[j]=1.0;
  }
  for(i=istart; i<iend; i++) {
    D->rtr[i]=1.0;
    D->rtd[i]=1.0;
    D->ltr[i]=1.0;
    D->ltd[i]=1.0;
  }

  if(D->pml==ON && D->pmlStart<iteration)  absorb_U(D);  else    ;


	j=jstart; m=0;
		for(i=istart; i<iend; i++) {
			D->BzNowR[m][i][j]=BzROld=D->BzR[m][i][j];
			D->EzNowR[m][i][j]=EzROld=D->EzR[m][i][j];
	
			D->BzR[m][i][j]=BzROld-dtBydr*(D->SrR[m][i-1][j]+D->SrR[m][i][j]+D->SlR[m][i-1][j]+D->SlR[m][i][j]);
			D->EzR[m][i][j]=EzROld+dtBydr*(D->PrR[m][i-1][j]+D->PrR[m][i][j]-D->PlR[m][i-1][j]-D->PlR[m][i][j])
								-M_PI*dt*(D->JzR[m][i-1][j]+D->JzR[m][i][j])
								+0.5*dF*dtBydz*(D->FR[m][i+1][j]-D->FR[m][i-1][j]);
		}

	for(m=0; m<numMode; m++) 
		for(j=jstart+1; j<jend; j++) {
			upr=D->upr[j]; upd=D->upd[j];

			r=(double)(j-jstart);
			for(i=istart; i<iend; i++) {
		    upPrR=0.5*(D->PrR[m][i-1][j]  +D->PrR[m][i][j]);
		    dnPrR=0.5*(D->PrR[m][i-1][j-1]+D->PrR[m][i][j-1]);
		    upPrI=0.5*(D->PrI[m][i-1][j]  +D->PrI[m][i][j]);
				dnPrI=0.5*(D->PrI[m][i-1][j-1]+D->PrI[m][i][j-1]);
				upSrR=0.5*(D->SrR[m][i-1][j]  +D->SrR[m][i][j]);
				dnSrR=0.5*(D->SrR[m][i-1][j-1]+D->SrR[m][i][j-1]);
				upSrI=0.5*(D->SrI[m][i-1][j]  +D->SrI[m][i][j]);
		    dnSrI=0.5*(D->SrI[m][i-1][j-1]+D->SrI[m][i][j-1]);
				upPlR=0.5*(D->PlR[m][i-1][j]  +D->PlR[m][i][j]);
				dnPlR=0.5*(D->PlR[m][i-1][j-1]+D->PlR[m][i][j-1]);
				upPlI=0.5*(D->PlI[m][i-1][j]  +D->PlI[m][i][j]);
				dnPlI=0.5*(D->PlI[m][i-1][j-1]+D->PlI[m][i][j-1]);
				upSlR=0.5*(D->SlR[m][i-1][j]  +D->SlR[m][i][j]);
				dnSlR=0.5*(D->SlR[m][i-1][j-1]+D->SlR[m][i][j-1]);
				upSlI=0.5*(D->SlI[m][i-1][j]  +D->SlI[m][i][j]);
				dnSlI=0.5*(D->SlI[m][i-1][j-1]+D->SlI[m][i][j-1]);
				PrR=0.5*(upPrR+dnPrR); PlR=0.5*(upPlR+dnPlR); 
				PrI=0.5*(upPrI+dnPrI); PlI=0.5*(upPlI+dnPlI); 
				SrR=0.5*(upSrR+dnSrR); SlR=0.5*(upSlR+dnSlR); 
				SrI=0.5*(upSrI+dnSrI); SlI=0.5*(upSlI+dnSlI); 

			 D->BzNowR[m][i][j]=BzROld=D->BzR[m][i][j];
			 D->BzNowI[m][i][j]=BzIOld=D->BzI[m][i][j];
			 D->EzNowR[m][i][j]=EzROld=D->EzR[m][i][j];
			 D->EzNowI[m][i][j]=EzIOld=D->EzI[m][i][j];

		    tmpr=upr;     tmpd=upd;

		    tmp=0.5/r*m*dtBydr*(SlI-SrI)
					  +0.5/r*dtBydr*((r+0.5)*(upPrR-upPlR)-(r-0.5)*(dnPrR-dnPlR))
						-M_PI*dt*(D->JzR[m][i-1][j]+D->JzR[m][i][j])
						+0.5*dF*dtBydz*(D->FR[m][i+1][j]-D->FR[m][i-1][j]);
				D->EzR[m][i][j]=tmpd*(EzROld+tmp*tmpr);
				tmp=-0.5/r*m*dtBydr*(SlR-SrR)
						+0.5/r*dtBydr*((r+0.5)*(upPrI-upPlI)-(r-0.5)*(dnPrI-dnPlI))
						-M_PI*dt*(D->JzI[m][i-1][j]+D->JzI[m][i][j])
						+0.5*dF*dtBydz*(D->FI[m][i+1][j]-D->FI[m][i-1][j]);
				D->EzI[m][i][j]=tmpd*(EzIOld+tmp*tmpr);

				tmp=-0.5/r*m*dtBydr*(PlI+PrI)
						  -0.5/r*dtBydr*((r+0.5)*(upSrR+upSlR)-(r-0.5)*(dnSrR+dnSlR));
					 D->BzR[m][i][j]=tmpd*(BzROld+tmp*tmpr);
					 tmp=0.5/r*m*dtBydr*(PlR+PrR)
						  -0.5/r*dtBydr*((r+0.5)*(upSrI+upSlI)-(r-0.5)*(dnSrI+dnSlI));
					 D->BzI[m][i][j]=tmpd*(BzIOld+tmp*tmpr);

				}
		  }


	// //boundary for Ez
	// if(D->Period==OFF) {
	// 	if(myrank==0) {
	// 		i=istart;
	// 		for(m=0; m<numMode; m++) 
 	// 			for(j=jstart; j<jend; j++) {
 	// 				D->EzR[m][i][j]=0.0;
 	// 				D->EzI[m][i][j]=0.0;
	//  			}
	// 	} else ;
	// 	if(myrank==nTasks-1) {
	// 		i=iend-1;
	// 		for(m=0; m<numMode; m++) 
	// 			for(j=jstart; j<jend; j++) {
	// 				D->EzR[m][i][j]=0.0;
	// 				D->EzI[m][i][j]=0.0;
	// 			}
	// 	} else ;

}

void PS_solve_Split(Domain *D,int iteration)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,dz,dt,dF;
  double rr,rd,upL,LdU,tmpr,tmpd,tmp,upr,upd;

  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr; dz=D->dz; dt=D->dt; dF=D->dF;
  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  rr=D->pmlr;    rd=D->pmld;

  for(j=jstart; j<jend; j++) {
    D->upr[j]=1.0;
    D->upd[j]=1.0;
  }
  for(i=istart; i<iend; i++) {
    D->rtr[i]=1.0;
    D->rtd[i]=1.0;
    D->ltr[i]=1.0;
    D->ltd[i]=1.0;
  }

  if(D->pml==ON && D->pmlStart<iteration)  absorb_U(D);  else    ;


	 for(m=0; m<numMode; m++) { 
		  for(j=jstart; j<jend; j++) {
				upr=D->upr[j]; upd=D->upd[j];
				r=(double)(j-jstart);
				tmpr=upr;     tmpd=upd;

				for(i=istart; i<iend; i++) {
					tmp=-0.5*dtBydr/(r+0.5)*m*(D->BzI[m][i+1][j+1]+D->BzI[m][i+1][j])
					  -dtBydr*(D->EzR[m][i+1][j+1]-D->EzR[m][i+1][j])
					  -2.0*M_PI*dt*D->JrR[m][i+1][j]
					  +dF*dtBydr*(D->FR[m][i+1][j+1]-D->FR[m][i+1][j]);
		      D->PlR[m][i][j]=tmpd*(D->PlR[m][i+1][j]+tmp*tmpr);
					tmp=0.5*dtBydr/(r+0.5)*m*(D->BzR[m][i+1][j+1]+D->BzR[m][i+1][j])
					  -dtBydr*(D->EzI[m][i+1][j+1]-D->EzI[m][i+1][j])
					  -2.0*M_PI*dt*D->JrI[m][i+1][j]
					  +dF*dtBydr*(D->FI[m][i+1][j+1]-D->FI[m][i+1][j]);
					D->PlI[m][i][j]=tmpd*(D->PlI[m][i+1][j]+tmp*tmpr);
					tmp=0.5*dtBydr/(r+0.5)*m*(D->EzI[m][i+1][j+1]+D->EzI[m][i+1][j])
					  -dtBydr*(D->BzR[m][i+1][j+1]-D->BzR[m][i+1][j])
					  -M_PI*dt*(D->JpR[m][i+1][j]+D->JpR[m][i+1][j+1])
					  -0.5*dF*m/(r+0.5)*dtBydr*(D->FI[m][i+1][j+1]+D->FI[m][i+1][j]);
					D->SlR[m][i][j]=tmpd*(D->SlR[m][i+1][j]+tmp*tmpr);
					tmp=-0.5*dtBydr/(r+0.5)*m*(D->EzR[m][i+1][j+1]+D->EzR[m][i+1][j])
					  -dtBydr*(D->BzI[m][i+1][j+1]-D->BzI[m][i+1][j])
					  -M_PI*dt*(D->JpI[m][i+1][j]+D->JpI[m][i+1][j+1])
					  +0.5*dF*m/(r+0.5)*dtBydr*(D->FR[m][i+1][j+1]+D->FR[m][i+1][j]);
					D->SlI[m][i][j]=tmpd*(D->SlI[m][i+1][j]+tmp*tmpr);
				}

				for(i=iend-1; i>=istart; i--) {
					 tmp=-0.5*dtBydr/(r+0.5)*m*(D->BzI[m][i][j]+D->BzI[m][i][j+1])
						  +dtBydr*(D->EzR[m][i][j+1]-D->EzR[m][i][j])
						  -2.0*M_PI*dt*D->JrR[m][i][j]
						  +dF*dtBydr*(D->FR[m][i][j+1]-D->FR[m][i][j]);
					 D->PrR[m][i][j]=tmpd*(D->PrR[m][i-1][j]+tmp*tmpr);
					 tmp=0.5*dtBydr/(r+0.5)*m*(D->BzR[m][i][j]+D->BzR[m][i][j+1])
						  +dtBydr*(D->EzI[m][i][j+1]-D->EzI[m][i][j])
						  -2.0*M_PI*dt*D->JrI[m][i][j]
						  +dF*dtBydr*(D->FI[m][i][j+1]-D->FI[m][i][j]);
					 D->PrI[m][i][j]=tmpd*(D->PrI[m][i-1][j]+tmp*tmpr);

					 tmp=-0.5*dtBydr/(r+0.5)*m*(D->EzI[m][i][j]+D->EzI[m][i][j+1])
						  -dtBydr*(D->BzR[m][i][j+1]-D->BzR[m][i][j])
						  -M_PI*dt*(D->JpR[m][i][j]+D->JpR[m][i][j+1])
						  -0.5*dF*m/(r+0.5)*dtBydr*(D->FI[m][i][j+1]+D->FI[m][i][j]);
					 D->SrR[m][i][j]=tmpd*(D->SrR[m][i-1][j]+tmp*tmpr);
					 tmp=0.5*dtBydr/(r+0.5)*m*(D->EzR[m][i][j]+D->EzR[m][i][j+1])
						  -dtBydr*(D->BzI[m][i][j+1]-D->BzI[m][i][j])
						  -M_PI*dt*(D->JpI[m][i][j]+D->JpI[m][i][j+1])
						  +0.5*dF*m/(r+0.5)*dtBydr*(D->FR[m][i][j+1]+D->FR[m][i][j]);
					 D->SrI[m][i][j]=tmpd*(D->SrI[m][i-1][j]+tmp*tmpr);
				}       //End of i

		  }         //End of j
	 }

}

void EzBz_Now_Split(Domain *D,int iteration)
{
  int i,j,m,numMode,istart,iend,jstart,jend;
  double BzROld,BzIOld,EzROld,EzIOld;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  for(m=0; m<numMode; m++)
     for(i=0; i<iend+3; i++)
        for(j=0; j<jend+3; j++) {
           D->EzNowR[m][i][j]=0.5*(D->EzR[m][i][j]+D->EzNowR[m][i][j]);
           D->EzNowI[m][i][j]=0.5*(D->EzI[m][i][j]+D->EzNowI[m][i][j]);
           D->BzNowR[m][i][j]=0.5*(D->BzR[m][i][j]+D->BzNowR[m][i][j]);
           D->BzNowI[m][i][j]=0.5*(D->BzI[m][i][j]+D->BzNowI[m][i][j]);
        }
}

void Yee_Now(Domain *D,int iteration)
{
  int i,j,m,numMode,istart,iend,jstart,jend;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  for(m=0; m<numMode; m++)
     for(i=0; i<iend+3; i++)
        for(j=0; j<jend+3; j++) {
           D->BrNowR[m][i][j]=0.5*(D->BrR[m][i][j]+D->BrNowR[m][i][j]);
           D->BrNowI[m][i][j]=0.5*(D->BrI[m][i][j]+D->BrNowI[m][i][j]);
           D->BpNowR[m][i][j]=0.5*(D->BpR[m][i][j]+D->BpNowR[m][i][j]);
           D->BpNowI[m][i][j]=0.5*(D->BpI[m][i][j]+D->BpNowI[m][i][j]);
           D->BzNowR[m][i][j]=0.5*(D->BzR[m][i][j]+D->BzNowR[m][i][j]);
           D->BzNowI[m][i][j]=0.5*(D->BzI[m][i][j]+D->BzNowI[m][i][j]);
        }
}

