#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

double maximum(double x1,double x2);
double minimum(double x1,double x2);

void updateCurrent1D_Split_1st(Domain *D,int nSpecies);
void updateCurrent2D_Split_1st(Domain *D,int nSpecies,int iteration);
void updateCurrent2D_Split_2nd(Domain *D);
void updateCurrent2D_Split_3rd(Domain *D);
void updateCurrent3D_Split_1st(Domain *D);
void updateCurrent3D_Split_2nd(Domain *D);
void updateCurrent3D_Split_3rd(Domain *D);

void updateCurrent_Split(Domain D,int iteration)
{
  void removeBoostIon();

  int nSpecies;

  nSpecies=D.nSpecies;
  if(D.boostOn==ON && D.boostIon==OFF)
    nSpecies=1;    

  switch((D.currentType-1)*3+D.dimension)  {
  case ((FIRST-1)*3+1) :    
    updateCurrent1D_Split_1st(&D,nSpecies);
    break;
  
  // 2D
  case ((FIRST-1)*3+2) :    
    updateCurrent2D_Split_1st(&D,nSpecies,iteration);
    break;

  case ((SECOND-1)*3+2) :
    updateCurrent2D_Split_2nd(&D);
    break;

  case ((THIRD-1)*3+2) :
    updateCurrent2D_Split_3rd(&D);
    break;

  //3D
  case ((FIRST-1)*3+3) :
    updateCurrent3D_Split_1st(&D);  
    break;

  case ((2-1)*3+3) :
    updateCurrent3D_Split_2nd(&D); 
    break;

  case ((3-1)*3+3) :
    updateCurrent3D_Split_3rd(&D);
    break;

  default :
    printf("In updateCurrent, what currentType(%d)?, what dimension(%d)?\n",D.currentType,D.dimension);
  }

  D.shareF[0]=D.Jx;
  D.shareF[1]=D.Jy;
  D.shareF[2]=D.Jz;
  MPI_TransferJ_Xplus(&D,3);
  MPI_TransferJ_Xminus(&D,3);
  MPI_TransferJ_Yplus(&D,3);
  MPI_TransferJ_Yminus(&D,3);    
  MPI_TransferJ_Zplus(&D,3);
  MPI_TransferJ_Zminus(&D,3); 
	if(D.Period==ON) {
    MPI_TransferJ_Period_X(&D,3);
    MPI_TransferJ_Period_Y(&D,3);	 
    MPI_TransferJ_Period_Z(&D,3);	 
	} else ;  
  D.Jx=D.shareF[0];
  D.Jy=D.shareF[1];
  D.Jz=D.shareF[2];

  //current share  
  D.shareF[0]=D.Jx;
  D.shareF[1]=D.Jy;
  D.shareF[2]=D.Jz;
  MPI_TransferF_Xminus(&D,3);
	if(D.Period==ON) {
    MPI_TransferF_Period_X(&D,3);
	} else ;
  D.Jx=D.shareF[0];
  D.Jy=D.shareF[1];
  D.Jz=D.shareF[2];
}

void updateCurrent1D_Split_1st(Domain *D,int nSpecies)
{
    int i,j,k,s,i1,i2;
    int istart,iend;
    int nxSub;
    double inverDt,x1,x2,y1,y2,xr,yr,zr,z1,z2;
    double Fx1,Fx2,Wx1,Wx2,dx,dt,wTimesQ;
    double wx,wy,vy,vz,invGam,xcc,xc;
    int intXc;
    
    ptclList *p;
    LoadList *LL;   
    Particle ***particle;
    particle=D->particle;

//    double maximum();
//    double minimum();

    double coeff[nSpecies];

    istart=D->istart;    iend=D->iend;    nxSub=D->nxSub;
    dt=D->dt;    dx=D->dx;   inverDt=1.0/D->dt;

    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       coeff[s]=LL->density/LL->criticalDensity;
       LL=LL->next;
       s++;
    }
    int myrank,nTasks;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    j=k=0;
    //initialize J
    if(D->fieldIonizationONOFF==OFF) {
        for(i=0; i<nxSub+5; i++)        {
          D->Jx[i][j][k]=0.0;
          D->Jy[i][j][k]=0.0;
          D->Jz[i][j][k]=0.0;
        }
    } else ;

      for(i=istart; i<iend; i++)
          for(s=0; s<nSpecies; s++)
          {
            p=particle[i][j][k].head[s]->pt;     
            while(p) 
            {
              invGam=1.0/sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
              wTimesQ=p->weight*p->charge;
              x2=p->x+i;             
              x1=p->oldX;          
              i1=(int)x1;         
              i2=(int)x2;   
              vy=p->p2*invGam;
              vz=p->p3*invGam;
              xc=0.5*(x1+x2);  
              intXc=(int)xc;  
              xcc=xc-intXc;  

              if(i1==i2) 
                xr=0.5*(x1+x2);
              else 
                xr=maximum(i1*1.0,i2*1.0);

              Fx1=(xr-x1);
              Fx2=(x2-xr);
              Wx1=0.5*(x1+xr)-i1;
              Wx2=0.5*(xr+x2)-i2;
   
              D->Jx[i1][j][k]    +=Fx1*coeff[s]*wTimesQ;
 
              D->Jx[i2][j][k]    +=Fx2*coeff[s]*wTimesQ;

              wx=1.0-xcc; 
              D->Jy[intXc][j][k]    +=wx*vy*coeff[s]*wTimesQ;
              D->Jz[intXc][j][k]    +=wx*vz*coeff[s]*wTimesQ;
              wx=xcc; 
              D->Jy[intXc+1][j][k]+=wx*vy*coeff[s]*wTimesQ;
              D->Jz[intXc+1][j][k]+=wx*vz*coeff[s]*wTimesQ;

              p=p->next;
            }    //End of while(p)
          }  	 //End of for(s)     
}

void updateCurrent2D_Split_1st(Domain *D,int nSpecies,int iteration)
{
	int i,j,k,n,s,ii,jj,i1,i2,j1,j2,k1,k2,gridI,gridJ,dataI[2],dataJ[2];
	int istart,iend,jstart,jend,kstart,kend;
	int nxSub,nySub,nzSub;
	double x1,x2,y1,y2,xr,yr,zr,z1,z2,dxBydt,dyBydt,dtOverdx,dtOverdy;
	double Fx,Wx[2],Fy,Wy[2],dataX[3],dataY[3];
	double vz,xc,yc,wTimesQ,posX1,posX2,posY1,posY2,invGam;
    
	ptclList *p;
	LoadList *LL;   
	Particle ***particle;
	particle=D->particle;

	double coeff[nSpecies];

	istart=D->istart;    iend=D->iend;
	jstart=D->jstart;    jend=D->jend;
	kstart=D->kstart;    kend=D->kend;

	nxSub=D->nxSub;    nySub=D->nySub;    nzSub=D->nzSub; 
	dxBydt=D->dx/D->dt; 		dyBydt=D->dy/D->dt;
	dtOverdx=D->dt/D->dx;	dtOverdy=D->dt/D->dy;
  
	s=0;
	LL=D->loadList;
	while(LL->next) {
		if(LL->species==2) coeff[s]=0.0;
		else coeff[s]=LL->density/LL->criticalDensity;
		LL=LL->next;
		s++;
	}

	int myrank,nTasks;
	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	k=0;
	//initialize J
	if(D->fieldIonizationONOFF==OFF) {
		for(i=0; i<nxSub+5; i++)
			for(j=0; j<nySub+5; j++)   {
            D->Jx[i][j][k]=0.0;
            D->Jy[i][j][k]=0.0;
            D->Jz[i][j][k]=0.0;
          }
  	} else ;

  	for(i=istart; i<iend; i++)
    	for(j=jstart; j<jend; j++)
      	for(s=0; s<nSpecies; s++)
      	{
        		p=particle[i][j][k].head[s]->pt;     
        		while(p) 
        		{
        			invGam=1.0/sqrt(1+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);					  
        			wTimesQ=p->charge*p->weight;
        			vz=p->p3*invGam;

					//------------------- Jy, Jz -------------------				  
          		x1=p->oldX;      y1=p->oldY;     
          		x2=p->x+i;       y2=p->y+j;      
					i1=(int)x1;      j1=(int)y1;     
          		i2=(int)x2;      j2=(int)y2;
	
          		if(i1==i2)        xr=0.5*(x1+x2);
          		else              xr=maximum(i1*1.0,i2*1.0);
          		if(j1==j2)        yr=0.5*(y1+y2);
          		else              yr=maximum(j1*1.0,j2*1.0);
	
          		dataX[0]=x1;   dataY[0]=y1;
          		dataX[1]=xr;   dataY[1]=yr;
          		dataX[2]=x2;   dataY[2]=y2;               
          		dataI[0]=i1;   dataJ[0]=j1;
          		dataI[1]=i2;   dataJ[1]=j2;               

					for(n=0; n<2; n++) {
            		gridI=dataI[n];   gridJ=dataJ[n];
            		posX1=dataX[n];   posX2=dataX[n+1];
            		posY1=dataY[n];   posY2=dataY[n+1];                  
						xc=0.5*(posX1+posX2);  yc=0.5*(posY1+posY2);                              
              
            		Fy=(posY2-posY1)*dyBydt;
            		Wx[1]=xc-gridI;   Wx[0]=1-Wx[1];
            		Wy[1]=yc-gridJ;   Wy[0]=1-Wy[1];
            		for(ii=0; ii<2; ii++) 
            			D->Jy[gridI+ii][gridJ][k]+=Fy*Wx[ii]*coeff[s]*wTimesQ;

            		for(ii=0; ii<2; ii++)
            			for(jj=0; jj<2; jj++) 
            		   	D->Jz[gridI+ii][gridJ+jj][k]+=Wx[ii]*Wy[jj]*vz*coeff[s]*wTimesQ;

         		}

					//------------------- Jx -------------------				                       
         		x1=(p->oldX+i+p->x)*0.5;   y1=(p->oldY+j+p->y)*0.5; 
               x2=x1+p->p1*invGam*dtOverdx;
               y2=y1+p->p2*invGam*dtOverdy;
					i1=(int)x1;      j1=(int)y1;     
          		i2=(int)x2;      j2=(int)y2;
               
          		if(i1==i2)        xr=0.5*(x1+x2);
          		else              xr=maximum(i1*1.0,i2*1.0);
          		if(j1==j2)        yr=0.5*(y1+y2);
          		else              yr=maximum(j1*1.0,j2*1.0);

          		dataX[0]=x1;   dataY[0]=y1;
          		dataX[1]=xr;   dataY[1]=yr;
          		dataX[2]=x2;   dataY[2]=y2;               
          		dataI[0]=i1;   dataJ[0]=j1;
          		dataI[1]=i2;   dataJ[1]=j2;               
          
          		for(n=0; n<2; n++) {
          		  	gridI=dataI[n];   gridJ=dataJ[n];
          		  	posX1=dataX[n];   posX2=dataX[n+1];
          		  	posY1=dataY[n];   posY2=dataY[n+1];                  

          		  	xc=0.5*(posX1+posX2);  yc=0.5*(posY1+posY2);                              
          		  	Wy[1]=yc-gridJ;   Wy[0]=1-Wy[1];

          		  	Fx=(posX2-posX1)*dxBydt;
          		  	for(jj=0; jj<2; jj++)
          		   	D->Jx[gridI][gridJ+jj][k]+=Fx*Wy[jj]*coeff[s]*wTimesQ;
          		}

          	p=p->next;
        	}    //End of while(p)
      }  	 //End of for(s)     

}

void updateCurrent2D_Split_2nd(Domain *D)
{
  	int i,j,k,s,n,i1,i2,j1,j2,k1,k2,ii,jj,kk,gridI,gridJ,dataI[2],dataJ[2];
  	int istart,iend,jstart,jend,kstart,kend;
  	int nxSub,nySub,nzSub;
  	double dxBydt,dyBydt,dzBydt,xr,yr,zr,wTimesQ,dtOverdx,dtOverdy,invGam;
  	double Fx[2],Fy[2],Fz[2],Wx[3],Wy[3],Wz[3],dataX[3],dataY[3];
  	double x1,x2,y1,y2,x,y,z,xc,yc,vz,posX1,posX2,posY1,posY2;
  	ptclList *p;
  	LoadList *LL;   
  	Particle ***particle;
  	particle=D->particle;
	
  	double coeff[D->nSpecies];

  	istart=D->istart; iend=D->iend;
  	jstart=D->jstart; jend=D->jend;
  	kstart=D->kstart; kend=D->kend;

  	nxSub=D->nxSub; nySub=D->nySub; nzSub=D->nzSub;
  	dxBydt=D->dx/D->dt; dyBydt=D->dy/D->dt;
	dtOverdx=D->dt/D->dx;	dtOverdy=D->dt/D->dy;	  

  	s=0;
  	LL=D->loadList;
  	while(LL->next)	{
  	  	if(LL->species==2) coeff[s]=0.0;
  	  	else coeff[s]=LL->density/LL->criticalDensity;
  	  	LL=LL->next;
  	  	s++;
  	}

  	int myrank,nTasks;
  	MPI_Status status;

  	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  	//initialize J
  	k=0;
  	if(D->fieldIonizationONOFF==OFF) {
  	   for(i=0; i<nxSub+5; i++)
  	      for(j=0; j<nySub+5; j++)   {
  	      	D->Jx[i][j][k]=0.0;
  	      	D->Jy[i][j][k]=0.0;
  	      	D->Jz[i][j][k]=0.0;
  	      }
  	} else ;

  	for(i=istart; i<iend; i++)
    	for(j=jstart; j<jend; j++)
      	for(s=0; s<D->nSpecies; s++)
      	{
        		p=particle[i][j][k].head[s]->pt;     
      
        		while(p) 
        		{
        			invGam=1.0/sqrt(1+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);					  
        			wTimesQ=p->charge*p->weight;
        			vz=p->p3*invGam;					 

					//------------------- Jy -------------------				        
          		x1=p->oldX;   y1=p->oldY;
          		x2=p->x+i;    y2=p->y+j;
                    
          		i1=(int)(x1+0.5);   j1=(int)(y1+0.5);
          		i2=(int)(x2+0.5);   j2=(int)(y2+0.5);

          		if(i1==i2)             xr=0.5*(x1+x2);
         		else                   xr=(i1+i2)*0.5;
         		if(j1==j2)             yr=0.5*(y1+y2);
         		else                   yr=(j1+j2)*0.5;

         		dataX[0]=x1;   dataY[0]=y1;
         		dataX[1]=xr;   dataY[1]=yr;
         		dataX[2]=x2;   dataY[2]=y2;               
         		dataI[0]=i1;   dataJ[0]=j1;
         		dataI[1]=i2;   dataJ[1]=j2;                  

					for(n=0; n<2; n++) {
            		gridI=dataI[n];   gridJ=dataJ[n];
            		posX1=dataX[n];   posX2=dataX[n+1];
            		posY1=dataY[n];   posY2=dataY[n+1];                  
            		xc=0.5*(posX1+posX2);  yc=0.5*(posY1+posY2);                  
            		x=xc-gridI;            y=yc-gridJ;

            		Fy[0]=(posY2-posY1)*(0.5-y)*dyBydt;
            		Fy[1]=(posY2-posY1)*(0.5+y)*dyBydt;
            		Wx[0]=0.5*(0.5-x)*(0.5-x);
            		Wx[1]=0.75-x*x;
            		Wx[2]=0.5*(0.5+x)*(0.5+x);
            		for(ii=0; ii<3; ii++)
            			for(jj=0; jj<2; jj++)
            		   	D->Jy[gridI-1+ii][gridJ-1+jj][k]+=Fy[jj]*Wx[ii]*coeff[s]*wTimesQ;

            		Wy[0]=0.5*(y-0.5)*(y-0.5);
            		Wy[1]=0.75-y*y;
            		Wy[2]=0.5*(y+0.5)*(y+0.5);
            		for(ii=0; ii<3; ii++)
            			for(jj=0; jj<3; jj++)
            		   	D->Jz[gridI-1+ii][gridJ-1+jj][k]+=Wx[ii]*Wy[jj]*vz*coeff[s]*wTimesQ;
         		}

					//------------------- Jx -------------------				        
         		x1=(p->oldX+i+p->x)*0.5;   y1=(p->oldY+j+p->y)*0.5; 
               x2=x1+p->p1*invGam*dtOverdx;
               y2=y1+p->p2*invGam*dtOverdy;
                    
          		i1=(int)(x1+0.5);   j1=(int)(y1+0.5);
          		i2=(int)(x2+0.5);   j2=(int)(y2+0.5);

          		if(i1==i2)             xr=0.5*(x1+x2);
          		else                   xr=(i1+i2)*0.5;
          		if(j1==j2)             yr=0.5*(y1+y2);
          		else                   yr=(j1+j2)*0.5;
					
          		dataX[0]=x1;   dataY[0]=y1;
          		dataX[1]=xr;   dataY[1]=yr;
          		dataX[2]=x2;   dataY[2]=y2;               
          		dataI[0]=i1;   dataJ[0]=j1;
          		dataI[1]=i2;   dataJ[1]=j2;   

					for(n=0; n<2; n++) {
            		gridI=dataI[n];   gridJ=dataJ[n];
            		posX1=dataX[n];   posX2=dataX[n+1];
            		posY1=dataY[n];   posY2=dataY[n+1];                  
            		xc=0.5*(posX1+posX2);  yc=0.5*(posY1+posY2);                  
            		x=xc-gridI;            y=yc-gridJ;

            		Fx[0]=(posX2-posX1)*(0.5-x)*dxBydt;
            		Fx[1]=(posX2-posX1)*(0.5+x)*dxBydt;
            		Wy[0]=0.5*(0.5-y)*(0.5-y);
            		Wy[1]=0.75-y*y;
            		Wy[2]=0.5*(0.5+y)*(0.5+y);
            		for(ii=0; ii<2; ii++)
            		  	for(jj=0; jj<3; jj++)
            		    	D->Jx[gridI-1+ii][gridJ-1+jj][k]+=Fx[ii]*Wy[jj]*coeff[s]*wTimesQ;
         		}

          		p=p->next;
        		}	//End of while(p)
      	}		//End of for(s)     

}

void updateCurrent2D_Split_3rd(Domain *D)
{
  	int i,j,k,s,n,i1,i2,j1,j2,k1,k2,ii,jj,kk,gridI,gridJ,dataI[2],dataJ[2];
  	int istart,iend,jstart,jend,kstart,kend;
  	int nxSub,nySub,nzSub;
  	double dxBydt,dyBydt,xr,yr,zr,x,y,z,dt,vz,wTimesQ,dtOverdx,dtOverdy,invGam;
  	double Fx[3],Fy[3],Fz[3],Wx[4],Wy[4],Wz[4],dataX[3],dataY[3];;
  	double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,xc,yc,posX1,posX2,posY1,posY2;
  	ptclList *p;
  	LoadList *LL;   
  	Particle ***particle;
  	particle=D->particle;
  	double coeff[D->nSpecies];

  	istart=D->istart;    iend=D->iend;
  	jstart=D->jstart;    jend=D->jend;
  	kstart=D->kstart;    kend=D->kend;

  	nxSub=D->nxSub;  nySub=D->nySub;  nzSub=D->nzSub;  
  	dxBydt=D->dx/D->dt; dyBydt=D->dy/D->dt;
	dtOverdx=D->dt/D->dx;	dtOverdy=D->dt/D->dy;	  
	
  	s=0;
  	LL=D->loadList;
  	while(LL->next)
  	{
  	  	if(LL->species==2) coeff[s]=0.0;
  	  	else coeff[s]=LL->density/LL->criticalDensity;
  	  	LL=LL->next;
  	  	s++;
  	}

  	int myrank,nTasks;
  	MPI_Status status;

  	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  	//initialize J
  	k=0;
  	if(D->fieldIonizationONOFF==OFF) {
  	   for(i=0; i<nxSub+5; i++)
  	      for(j=0; j<nySub+5; j++)   {
  	      	D->Jx[i][j][k]=0.0;
  	      	D->Jy[i][j][k]=0.0;
  	      	D->Jz[i][j][k]=0.0;
  	      }
  	} else ;

  	for(i=istart; i<iend; i++)
    	for(j=jstart; j<jend; j++)
      	for(s=0; s<D->nSpecies; s++)
      	{
        		p=particle[i][j][k].head[s]->pt;     
      
        		while(p) 
        		{
        			invGam=1.0/sqrt(1+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);					  
        			wTimesQ=p->charge*p->weight;
        			vz=p->p3*invGam;		

					//------------------- Jy, Jz -------------------				        
          		x1=p->oldX;   y1=p->oldY;
          		x2=p->x+i;    y2=p->y+j;

          		i1=(int)(x1);   j1=(int)(y1);
          		i2=(int)(x2);   j2=(int)(y2);
					
          		if(i1==i2)             xr=0.5*(x1+x2);
          		else                   xr=maximum(i1*1.0,i2*1.0);
          		if(j1==j2)             yr=0.5*(y1+y2);
          		else                   yr=maximum(j1*1.0,j2*1.0);

          		dataX[0]=x1;   dataY[0]=y1;
          		dataX[1]=xr;   dataY[1]=yr;
          		dataX[2]=x2;   dataY[2]=y2;               
          		dataI[0]=i1;   dataJ[0]=j1;
          		dataI[1]=i2;   dataJ[1]=j2;                  

					for(n=0; n<2; n++) {
            		gridI=dataI[n];   gridJ=dataJ[n];
            		posX1=dataX[n];   posX2=dataX[n+1];
            		posY1=dataY[n];   posY2=dataY[n+1];                  

            		xc=0.5*(posX1+posX2);  yc=0.5*(posY1+posY2);                  
            		x=xc-gridI;            y=yc-gridJ;

            		Fy[0]=(posY2-posY1)*0.5*(1-y)*(1-y)*dyBydt;
            		Fy[1]=(posY2-posY1)*(0.75-(0.5-y)*(0.5-y))*dyBydt;
            		Fy[2]=(posY2-posY1)*0.5*y*y*dyBydt;
            		x1=1+x;
            		x2=x;
            		x3=1-x;
            		x4=2-x;
            		Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
            		Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
            		Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
            		Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
            		for(ii=0; ii<4; ii++)
            		  for(jj=0; jj<3; jj++)
            		    D->Jy[gridI-1+ii][gridJ-1+jj][k]+=Fy[jj]*Wx[ii]*coeff[s]*wTimesQ;

            		x1=1+x;
            		x2=x;
            		x3=1-x;
            		x4=2-x;
            		Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
            		Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
            		Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
            		Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
            		y1=1+y;
            		y2=y;
            		y3=1-y;
            		y4=2-y;
            		Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
            		Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
            		Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
            		Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;

            		for(ii=0; ii<4; ii++)
            		  for(jj=0; jj<4; jj++)
            		    D->Jz[gridI-1+ii][gridJ-1+jj][k]+=Wx[ii]*Wy[jj]*vz*coeff[s]*wTimesQ;
         		}

					//------------------- Jx -------------------				        
         		x1=(p->oldX+i+p->x)*0.5;   y1=(p->oldY+j+p->y)*0.5; 
               x2=x1+p->p1*invGam*dtOverdx;
               y2=y1+p->p2*invGam*dtOverdy;

					i1=(int)(x1);   j1=(int)(y1);
					i2=(int)(x2);   j2=(int)(y2);

					if(i1==i2)             xr=0.5*(x1+x2);
					else                   xr=maximum(i1*1.0,i2*1.0);
					if(j1==j2)             yr=0.5*(y1+y2);
					else                   yr=maximum(j1*1.0,j2*1.0);					

					dataX[0]=x1;   dataY[0]=y1;
					dataX[1]=xr;   dataY[1]=yr;
					dataX[2]=x2;   dataY[2]=y2;               
					dataI[0]=i1;   dataJ[0]=j1;
					dataI[1]=i2;   dataJ[1]=j2;                  

					for(n=0; n<2; n++) {
            		gridI=dataI[n];   gridJ=dataJ[n];
            		posX1=dataX[n];   posX2=dataX[n+1];
            		posY1=dataY[n];   posY2=dataY[n+1];                  

            		xc=0.5*(posX1+posX2);  yc=0.5*(posY1+posY2);                  
            		x=xc-gridI;            y=yc-gridJ;

            		Fx[0]=(posX2-posX1)*0.5*(1-x)*(1-x)*dxBydt;
            		Fx[1]=(posX2-posX1)*(0.75-(0.5-x)*(0.5-x))*dxBydt;
            		Fx[2]=(posX2-posX1)*0.5*x*x*dxBydt;
            		y1=1+y;
            		y2=y;
            		y3=1-y;
            		y4=2-y;
            		Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
            		Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
            		Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
            		Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;

            		for(ii=0; ii<3; ii++)
            		  	for(jj=0; jj<4; jj++)
            		    	D->Jx[gridI-1+ii][gridJ-1+jj][k]+=Fx[ii]*Wy[jj]*coeff[s]*wTimesQ;
         		}
        			p=p->next;
      		}	//End of while(p)
    		}	//End of for(s)     

}


void updateCurrent3D_Split_1st(Domain *D)
{
   int i,j,k,s,n,ii,jj,kk,i1,i2,j1,j2,k1,k2,dataI[2],dataJ[2],dataK[2];
   int istart,iend,jstart,jend,kstart,kend,gridI,gridJ,gridK;
   double x1,x2,y1,y2,xr,yr,zr,z1,z2,wTimesQ,posX1,posX2,posY1,posY2,posZ1,posZ2;
   double Fx,Wx[2],Fy,Wy[2],Fz,Wz[2],dataX[3],dataY[3],dataZ[3];
   double dxBydt,dyBydt,dzBydt,xc,yc,zc,dtOverdx,dtOverdy,dtOverdz,invGam;
   ptclList *p;
   LoadList *LL;   
   Particle ***particle;
   particle=D->particle;

   double coeff[D->nSpecies];

   istart=D->istart;    iend=D->iend;
   jstart=D->jstart;    jend=D->jend;
   kstart=D->kstart;    kend=D->kend;

   dxBydt=D->dx/D->dt; 		dyBydt=D->dy/D->dt; 		dzBydt=D->dz/D->dt;
	dtOverdx=D->dt/D->dx;	dtOverdy=D->dt/D->dy;	dtOverdz=D->dt/D->dz;

   s=0;
   LL=D->loadList;
   while(LL->next)  {
     	if(LL->species==2) coeff[s]=0.0;
     	else coeff[s]=LL->density/LL->criticalDensity;
     	LL=LL->next;
     	s++;
   }

   int myrank,nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   //initialize J
   if(D->fieldIonizationONOFF==OFF) {
      for(i=0; i<iend+3; i++)
         for(j=0; j<jend+3; j++)
         	for(k=0; k<kend+3; k++)  {
            	D->Jx[i][j][k]=0.0;
            	D->Jy[i][j][k]=0.0;
            	D->Jz[i][j][k]=0.0;
         	}
   } else ;

   for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
         for(k=kstart; k<kend; k++)
            for(s=0; s<D->nSpecies; s++)
            {
               p=particle[i][j][k].head[s]->pt;     
         
               while(p) 
               {
                  wTimesQ=p->weight*p->charge;
						invGam=1.0/sqrt(1+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);

						//------------------- Jy, Jz ------------------- 
                  x1=p->oldX;      y1=p->oldY;      z1=p->oldZ;
                  x2=p->x+i;       y2=p->y+j;       z2=p->z+k;

                  i1=(int)x1;      j1=(int)y1;      k1=(int)z1;
                  i2=(int)x2;      j2=(int)y2;      k2=(int)z2;
              
                  if(i1==i2)          xr=0.5*(x1+x2);
                  else                xr=maximum(i1*1.0,i2*1.0);
                  if(j1==j2)          yr=0.5*(y1+y2);
                  else                yr=maximum(j1*1.0,j2*1.0);
                  if(k1==k2)          zr=0.5*(z1+z2);
                  else                zr=maximum(k1*1.0,k2*1.0);

                  dataX[0]=x1;   dataX[1]=xr;   dataX[2]=x2;
                  dataY[0]=y1;   dataY[1]=yr;   dataY[2]=y2;
                  dataZ[0]=z1;   dataZ[1]=zr;   dataZ[2]=z2;                                    
                  dataI[0]=i1;   dataJ[0]=j1;   dataK[0]=k1;
                  dataI[1]=i2;   dataJ[1]=j2;   dataK[1]=k2;
                       
                  for(n=0; n<2; n++) {
                     gridI=dataI[n];   gridJ=dataJ[n];   gridK=dataK[n];
                     posX1=dataX[n];   posY1=dataY[n];   posZ1=dataZ[n];
                     posX2=dataX[n+1]; posY2=dataY[n+1]; posZ2=dataZ[n+1];                     
                     xc=0.5*(posX1+posX2);  yc=0.5*(posY1+posY2);  zc=0.5*(posZ1+posZ2);
                     
                     Wx[1]=xc-gridI;   Wx[0]=1-Wx[1];
                     Wy[1]=yc-gridJ;   Wy[0]=1-Wy[1];                 
                     Wz[1]=zc-gridK;   Wz[0]=1-Wz[1];                                      
                     
                     Fy=(posY2-posY1)*dyBydt;
                     for(ii=0; ii<2; ii++)
                        for(kk=0; kk<2; kk++)                  
                           D->Jy[gridI+ii][gridJ][gridK+kk]+=Fy*Wx[ii]*Wz[kk]*coeff[s]*wTimesQ;

                     Fz=(posZ2-posZ1)*dzBydt;
                     for(ii=0; ii<2; ii++)
                        for(jj=0; jj<2; jj++)                  
                           D->Jz[gridI+ii][gridJ+jj][gridK]+=Fz*Wx[ii]*Wy[jj]*coeff[s]*wTimesQ;
                  }

						//------------------- Jx------------------- 
                  x1=(p->oldX+i+p->x)*0.5;   y1=(p->oldY+j+p->y)*0.5;   z1=(p->oldZ+k+p->z)*0.5;
						x2=x1+p->p1*invGam*dtOverdx;
               	y2=y1+p->p2*invGam*dtOverdy;
						z2=z1+p->p3*invGam*dtOverdz;

                  i1=(int)x1;      j1=(int)y1;      k1=(int)z1;
                  i2=(int)x2;      j2=(int)y2;      k2=(int)z2;
              
                  if(i1==i2)          xr=0.5*(x1+x2);
                  else                xr=maximum(i1*1.0,i2*1.0);
                  if(j1==j2)          yr=0.5*(y1+y2);
                  else                yr=maximum(j1*1.0,j2*1.0);
                  if(k1==k2)          zr=0.5*(z1+z2);
                  else                zr=maximum(k1*1.0,k2*1.0);

                  dataX[0]=x1;   dataX[1]=xr;   dataX[2]=x2;
                  dataY[0]=y1;   dataY[1]=yr;   dataY[2]=y2;
                  dataZ[0]=z1;   dataZ[1]=zr;   dataZ[2]=z2;                                    
                  dataI[0]=i1;   dataJ[0]=j1;   dataK[0]=k1;
                  dataI[1]=i2;   dataJ[1]=j2;   dataK[1]=k2;
                       
                  for(n=0; n<2; n++) {
                     gridI=dataI[n];   gridJ=dataJ[n];   gridK=dataK[n];
                     posX1=dataX[n];   posY1=dataY[n];   posZ1=dataZ[n];
                     posX2=dataX[n+1]; posY2=dataY[n+1]; posZ2=dataZ[n+1];                     
                     xc=0.5*(posX1+posX2);  yc=0.5*(posY1+posY2);  zc=0.5*(posZ1+posZ2);
                     
                     Wx[1]=xc-gridI;   Wx[0]=1-Wx[1];
                     Wy[1]=yc-gridJ;   Wy[0]=1-Wy[1];                 
                     Wz[1]=zc-gridK;   Wz[0]=1-Wz[1];                                      
                     
                     Fx=(posX2-posX1)*dxBydt;
                     for(jj=0; jj<2; jj++)
                        for(kk=0; kk<2; kk++)                  
                           D->Jx[gridI][gridJ+jj][gridK+kk]+=Fx*Wy[jj]*Wz[kk]*coeff[s]*wTimesQ;
                  }
                
               p=p->next;
            }	//End of while(p)
         }	   //End of for(s)     
      
}


void updateCurrent3D_Split_2nd(Domain *D)
{
  	int i,j,k,s,n,ii,jj,kk,i1,i2,j1,j2,k1,k2,dataI[2],dataJ[2],dataK[2];
  	int istart,iend,jstart,jend,kstart,kend,gridI,gridJ,gridK;
  	double x1,x2,y1,y2,z1,z2,xr,yr,zr,x,y,z,wTimesQ,posX1,posX2,posY1,posY2,posZ1,posZ2;
  	double Fx[2],Fy[2],Fz[2],Wx[3],Wy[3],Wz[3],dataX[3],dataY[3],dataZ[3];
  	double dxBydt,dyBydt,dzBydt,xc,yc,zc,dtOverdx,dtOverdy,dtOverdz,invGam;

  	ptclList *p;
  	LoadList *LL;   
  	Particle ***particle;
  	particle=D->particle;

  	double coeff[D->nSpecies];

  	istart=D->istart;    iend=D->iend;
  	jstart=D->jstart;    jend=D->jend;
  	kstart=D->kstart;    kend=D->kend;

  	dxBydt=D->dx/D->dt; dyBydt=D->dy/D->dt; dzBydt=D->dz/D->dt;
	dtOverdx=D->dt/D->dx;	dtOverdy=D->dt/D->dy;	dtOverdz=D->dt/D->dz;

  	s=0;
  	LL=D->loadList;
  	while(LL->next)    {
  	  	if(LL->species==2) coeff[s]=0.0;
  	  	else coeff[s]=LL->density/LL->criticalDensity;
  	  	LL=LL->next;
  	  	s++;
  	}

  	int myrank,nTasks;
  	MPI_Status status;

  	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  //initialize J
  if(D->fieldIonizationONOFF==OFF) {
      for(i=0; i<iend+3; i++)
        	for(j=0; j<jend+3; j++)
          	for(k=0; k<kend+3; k++)  {
            	D->Jx[i][j][k]=0.0;
            	D->Jy[i][j][k]=0.0;
            	D->Jz[i][j][k]=0.0;
          	}
  	} else ;

  	for(i=istart; i<iend; i++)
    	for(j=jstart; j<jend; j++)
      	for(k=kstart; k<kend; k++)
        		for(s=0; s<D->nSpecies; s++)
        		{
          		p=particle[i][j][k].head[s]->pt;     

          		while(p) 
          		{
                  wTimesQ=p->weight*p->charge;
						invGam=1.0/sqrt(1+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);

						//------------------- Jy -------------------				        
            		x1=p->oldX;   y1=p->oldY;   z1=p->oldZ;
            		x2=p->x+i;    y2=p->y+j;    z2=p->z+k;

            		i1=(int)(x1+0.5);   j1=(int)(y1+0.5); k1=(int)(z1+0.5);
            		i2=(int)(x2+0.5);   j2=(int)(y2+0.5); k2=(int)(z2+0.5);

            		if(i1==i2)             xr=0.5*(x1+x2);
            		else                   xr=(i1+i2)*0.5;
            		if(j1==j2)             yr=0.5*(y1+y2);
            		else                   yr=(j1+j2)*0.5;
            		if(k1==k2)             zr=0.5*(z1+z2);
            		else                   zr=(k1+k2)*0.5;
		
            		dataX[0]=x1;   dataY[0]=y1; dataZ[0]=z1;
            		dataX[1]=xr;   dataY[1]=yr; dataZ[1]=zr;
            		dataX[2]=x2;   dataY[2]=y2; dataZ[2]=z2;               
            		dataI[0]=i1;   dataJ[0]=j1; dataK[0]=k1;
            		dataI[1]=i2;   dataJ[1]=j2; dataK[1]=k2;                  

			  			for(n=0; n<2; n++) {
              			gridI=dataI[n];   gridJ=dataJ[n];      gridK=dataK[n];   
              			posX1=dataX[n];   posX2=dataX[n+1];
              			posY1=dataY[n];   posY2=dataY[n+1];
              			posZ1=dataZ[n];   posZ2=dataZ[n+1];

              			xc=0.5*(posX1+posX2);  yc=0.5*(posY1+posY2);    zc=0.5*(posZ1+posZ2);                  
              			x=xc-gridI;            y=yc-gridJ;              z=zc-gridK;

              			Fy[0]=(posY2-posY1)*(0.5-y)*dyBydt;
              			Fy[1]=(posY2-posY1)*(0.5+y)*dyBydt;
              			Fz[0]=(posZ2-posZ1)*(0.5-z)*dzBydt;
              			Fz[1]=(posZ2-posZ1)*(0.5+z)*dzBydt;
              			Wx[0]=0.5*(0.5-x)*(0.5-x);
              			Wx[1]=0.75-x*x;
              			Wx[2]=0.5*(0.5+x)*(0.5+x);
              			Wy[0]=0.5*(0.5-y)*(0.5-y);
              			Wy[1]=0.75-y*y;
              			Wy[2]=0.5*(0.5+y)*(0.5+y);
              			Wz[0]=0.5*(0.5-z)*(0.5-z);
              			Wz[1]=0.75-z*z;
              			Wz[2]=0.5*(0.5+z)*(0.5+z);
              			for(ii=0; ii<3; ii++)
              			  	for(jj=0; jj<2; jj++)
              			    	for(kk=0; kk<3; kk++)
              			      	D->Jy[gridI-1+ii][gridJ-1+jj][gridK-1+kk]+=Fy[jj]*Wx[jj]*Wz[kk]*coeff[s]*wTimesQ;
              			for(ii=0; ii<3; ii++)
              			  	for(jj=0; jj<3; jj++)
              			    	for(kk=0; kk<2; kk++)
              			      	D->Jz[gridI-1+ii][gridJ-1+jj][gridK-1+kk]+=Fz[kk]*Wx[ii]*Wy[jj]*coeff[s]*wTimesQ;
            		}

						//------------------- Jx -------------------				        
                  x1=(p->oldX+i+p->x)*0.5;   y1=(p->oldY+j+p->y)*0.5;   z1=(p->oldZ+k+p->z)*0.5;
						x2=x1+p->p1*invGam*dtOverdx;
               	y2=y1+p->p2*invGam*dtOverdy;
						z2=z1+p->p3*invGam*dtOverdz;

            		i1=(int)(x1+0.5);   j1=(int)(y1+0.5); k1=(int)(z1+0.5);
            		i2=(int)(x2+0.5);   j2=(int)(y2+0.5); k2=(int)(z2+0.5);

            		if(i1==i2)             xr=0.5*(x1+x2);
            		else                   xr=(i1+i2)*0.5;
            		if(j1==j2)             yr=0.5*(y1+y2);
            		else                   yr=(j1+j2)*0.5;
            		if(k1==k2)             zr=0.5*(z1+z2);
            		else                   zr=(k1+k2)*0.5;

            		dataX[0]=x1;   dataY[0]=y1; dataZ[0]=z1;
            		dataX[1]=xr;   dataY[1]=yr; dataZ[1]=zr;
            		dataX[2]=x2;   dataY[2]=y2; dataZ[2]=z2;               
            		dataI[0]=i1;   dataJ[0]=j1; dataK[0]=k1;
            		dataI[1]=i2;   dataJ[1]=j2; dataK[1]=k2;                  

			  			for(n=0; n<2; n++) {
              			gridI=dataI[n];   gridJ=dataJ[n];      gridK=dataK[n];   
              			posX1=dataX[n];   posX2=dataX[n+1];
              			posY1=dataY[n];   posY2=dataY[n+1];
              			posZ1=dataZ[n];   posZ2=dataZ[n+1];

              			xc=0.5*(posX1+posX2);  yc=0.5*(posY1+posY2);    zc=0.5*(posZ1+posZ2);                  
              			x=xc-gridI;            y=yc-gridJ;              z=zc-gridK;                

              			Fx[0]=(posX2-posX1)*(0.5-x)*dxBydt;
              			Fx[1]=(posX2-posX1)*(0.5+x)*dxBydt;
              			Wy[0]=0.5*(0.5-y)*(0.5-y);
              			Wy[1]=0.75-y*y;
              			Wy[2]=0.5*(0.5+y)*(0.5+y);
              			Wz[0]=0.5*(0.5-z)*(0.5-z);
              			Wz[1]=0.75-z*z;
              			Wz[2]=0.5*(0.5+z)*(0.5+z);          

              			for(ii=0; ii<2; ii++)
              			  	for(jj=0; jj<3; jj++)
              			    	for(kk=0; kk<3; kk++)
                    				D->Jx[gridI-1+ii][gridJ-1+jj][gridK-1+kk]+=Fx[ii]*Wy[jj]*Wz[kk]*coeff[s]*wTimesQ;
            		}

            		p=p->next;
          		}	//End of while(p)
        		}		//End of for(s)     

}



void updateCurrent3D_Split_3rd(Domain *D)
{
  	int i,j,k,s,n,ii,jj,kk,i1,i2,j1,j2,k1,k2,dataI[2],dataJ[2],dataK[2];
  	int istart,iend,jstart,jend,kstart,kend,gridI,gridJ,gridK;
  	double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,xr,yr,zr,x,y,z,wTimesQ,posX1,posX2,posY1,posY2,posZ1,posZ2;
  	double Fx[3],Fy[3],Fz[3],Wx[4],Wy[4],Wz[4],dataX[3],dataY[3],dataZ[3];
  	double dxBydt,dyBydt,dzBydt,xc,yc,zc,dtOverdx,dtOverdy,dtOverdz,invGam;
  	ptclList *p;
  	LoadList *LL;   
  	Particle ***particle;
  	particle=D->particle;

  	double coeff[D->nSpecies];

  	istart=D->istart; iend=D->iend;
  	jstart=D->jstart; jend=D->jend;
  	kstart=D->kstart; kend=D->kend;

  	dxBydt=D->dx/D->dt; dyBydt=D->dy/D->dt; dzBydt=D->dz/D->dt;
	dtOverdx=D->dt/D->dx;	dtOverdy=D->dt/D->dy;	dtOverdz=D->dt/D->dz;  

  	s=0;
  	LL=D->loadList;
  	while(LL->next) {
  	  	if(LL->species==2) coeff[s]=0.0;
  	  	else coeff[s]=LL->density/LL->criticalDensity;
  	  	LL=LL->next;
  	  	s++;
  	}

  	int myrank,nTasks;
  	MPI_Status status;

  	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  	//initialize J
  	if(D->fieldIonizationONOFF==OFF) {
  	   for(i=0; i<iend+3; i++)
  	      for(j=0; j<jend+3; j++)
  	        	for(k=0; k<kend+3; k++)  {
  	         	D->Jx[i][j][k]=0.0;
  	         	D->Jy[i][j][k]=0.0;
  	         	D->Jz[i][j][k]=0.0;
  	        	}
  	} else ;

  	for(i=istart; i<iend; i++)
    	for(j=jstart; j<jend; j++)
      	for(k=kstart; k<kend; k++)
        		for(s=0; s<D->nSpecies; s++)
        		{
          		p=particle[i][j][k].head[s]->pt;     
       
          		while(p) 
          		{
                  wTimesQ=p->weight*p->charge;
						invGam=1.0/sqrt(1+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);

    					//------------------- Jy, Jz -------------------				        
            		x1=p->oldX;   y1=p->oldY;   z1=p->oldZ;
            		x2=p->x+i;    y2=p->y+j;    z2=p->z+k;

            		i1=(int)(x1);   j1=(int)(y1);   k1=(int)(z1);
            		i2=(int)(x2);   j2=(int)(y2);   k2=(int)(z2);

            		if(i1==i2)             xr=0.5*(x1+x2);
            		else                   xr=maximum(i1*1.0,i2*1.0);
            		if(j1==j2)             yr=0.5*(y1+y2);
            		else                   yr=maximum(j1*1.0,j2*1.0);
            		if(k1==k2)             zr=0.5*(z1+z2);
            		else                   zr=maximum(k1*1.0,k2*1.0);

            		dataX[0]=x1;   dataY[0]=y1;   dataZ[0]=z1;
            		dataX[1]=xr;   dataY[1]=yr;   dataZ[1]=zr;
            		dataX[2]=x2;   dataY[2]=y2;   dataZ[2]=z2;               
            		dataI[0]=i1;   dataJ[0]=j1;   dataK[0]=k1;
            		dataI[1]=i2;   dataJ[1]=j2;   dataK[1]=k2;              

				    	for(n=0; n<2; n++) {
              			gridI=dataI[n];   gridJ=dataJ[n];     gridK=dataK[n];
              			posX1=dataX[n];   posX2=dataX[n+1];
              			posY1=dataY[n];   posY2=dataY[n+1];                  
              			posZ1=dataZ[n];   posZ2=dataZ[n+1];  

              			xc=0.5*(posX1+posX2);  yc=0.5*(posY1+posY2);  zc=0.5*(posZ1+posZ2);               
              			x=xc-gridI;            y=yc-gridJ;            z=zc-gridK;

              			x1=1+x;
              			x2=x;
              			x3=1-x;
              			x4=2-x;
              			Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
              			Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
              			Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
              			Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
              			y1=1+y;
              			y2=y;
              			y3=1-y;
              			y4=2-y;
              			Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
              			Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
              			Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
              			Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;
              			z1=1+z;
              			z2=z;
              			z3=1-z;
              			z4=2-z;
              			Wz[0]=(2-z1)*(2-z1)*(2-z1)/6.0;
              			Wz[1]=(4-6*z2*z2+3*z2*z2*z2)/6.0;
              			Wz[2]=(4-6*z3*z3+3*z3*z3*z3)/6.0;
              			Wz[3]=(2-z4)*(2-z4)*(2-z4)/6.0;

              			Fy[0]=(posY2-posY1)*0.5*(1-y)*(1-y)*dyBydt;
              			Fy[1]=(posY2-posY1)*(0.75-(0.5-y)*(0.5-y))*dyBydt;
              			Fy[2]=(posY2-posY1)*0.5*y*y*dyBydt;
              			for(ii=0; ii<4; ii++)
              			  for(jj=0; jj<3; jj++)
              			    for(kk=0; kk<4; kk++)
              			      D->Jy[gridI-1+ii][gridJ-1+jj][gridK-1+kk]+=Fy[jj]*Wx[ii]*Wz[kk]*coeff[s]*wTimesQ;

              			Fz[0]=(posZ2-posZ1)*0.5*(1-z)*(1-z)*dzBydt;
              			Fz[1]=(posZ2-posZ1)*(0.75-(0.5-z)*(0.5-z))*dzBydt;
              			Fz[2]=(posZ2-posZ1)*0.5*z*z*dzBydt;   
              			for(ii=0; ii<4; ii++)
              			  for(jj=0; jj<4; jj++)
              			    for(kk=0; kk<3; kk++)
              			      D->Jz[gridI-1+ii][gridJ-1+jj][gridK-1+kk]+=Fz[kk]*Wx[ii]*Wy[jj]*coeff[s]*wTimesQ;
            		}

						//------------------- Jx -------------------	
                  x1=(p->oldX+i+p->x)*0.5;   y1=(p->oldY+j+p->y)*0.5;   z1=(p->oldZ+k+p->z)*0.5;
						x2=x1+p->p1*invGam*dtOverdx;
               	y2=y1+p->p2*invGam*dtOverdy;
						z2=z1+p->p3*invGam*dtOverdz;
                  
            		i1=(int)(x1);   j1=(int)(y1);   k1=(int)(z1);
            		i2=(int)(x2);   j2=(int)(y2);   k2=(int)(z2);

            		if(i1==i2)             xr=0.5*(x1+x2);
            		else                   xr=maximum(i1*1.0,i2*1.0);
            		if(j1==j2)             yr=0.5*(y1+y2);
            		else                   yr=maximum(j1*1.0,j2*1.0);
            		if(k1==k2)             zr=0.5*(z1+z2);
            		else                   zr=maximum(k1*1.0,k2*1.0);

            		dataX[0]=x1;   dataY[0]=y1;   dataZ[0]=z1;
            		dataX[1]=xr;   dataY[1]=yr;   dataZ[1]=zr;
            		dataX[2]=x2;   dataY[2]=y2;   dataZ[2]=z2;               
            		dataI[0]=i1;   dataJ[0]=j1;   dataK[0]=k1;
            		dataI[1]=i2;   dataJ[1]=j2;   dataK[1]=k2;              

				   	for(n=0; n<2; n++) {
              			gridI=dataI[n];   gridJ=dataJ[n];     gridK=dataK[n];
              			posX1=dataX[n];   posX2=dataX[n+1];
              			posY1=dataY[n];   posY2=dataY[n+1];                  
              			posZ1=dataZ[n];   posZ2=dataZ[n+1];  

              			xc=0.5*(posX1+posX2);  yc=0.5*(posY1+posY2);  zc=0.5*(posZ1+posZ2);               
              			x=xc-gridI;            y=yc-gridJ;            z=zc-gridK;

              			y1=1+y;
              			y2=y;
              			y3=1-y;
              			y4=2-y;
              			Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
              			Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
              			Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
              			Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;
              			z1=1+z;
              			z2=z;
              			z3=1-z;
              			z4=2-z;
              			Wz[0]=(2-z1)*(2-z1)*(2-z1)/6.0;
              			Wz[1]=(4-6*z2*z2+3*z2*z2*z2)/6.0;
              			Wz[2]=(4-6*z3*z3+3*z3*z3*z3)/6.0;
              			Wz[3]=(2-z4)*(2-z4)*(2-z4)/6.0;

              			Fx[0]=(posX2-posX1)*0.5*(1-x)*(1-x)*dxBydt;
              			Fx[1]=(posX2-posX1)*(0.75-(0.5-x)*(0.5-x))*dxBydt;
              			Fx[2]=(posX2-posX1)*0.5*x*x*dxBydt;
              			for(ii=0; ii<3; ii++)
              			  	for(jj=0; jj<4; jj++)
              			    	for(kk=0; kk<4; kk++)
              			      	D->Jx[gridI-1+ii][gridJ-1+jj][gridK-1+kk]+=Fx[ii]*Wy[jj]*Wz[kk]*coeff[s]*wTimesQ;
            		}

            		p=p->next;
          		}	//End of while(p)
        		}	//End of for(s)     

}



