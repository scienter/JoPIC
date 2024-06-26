#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_randist.h>


double gaussian_dist(double sig);
void MPI_TransferBeamDen_Xplus(Domain *D,double **f1,int ny,int share);
void MPI_TransferBeamDen_Xminus(Domain *D,double **f1,int ny,int share);
void MPI_TransferBeamField_Xplus(Domain *D,double **f1,int ny,int share);
void MPI_TransferBeamField_Xminus(Domain *D,double **f1,int ny,int share);
void loadBeamPlasma(Domain *D,LoadList *LL,int s,int iteration);

void loadBeam(Domain *D,LoadList *LL,int s,int iteration)
{
  switch(LL->type)  {
  case Polygon:
    ;
    break;
  case Beam:
    if(iteration==LL->loadingStep) { 
      D->shiftStart=D->nx;
      loadBeamPlasma(D,LL,s,iteration);
    } else ;
    break;
  default:
    ;
  }
}


void loadBeamPlasma(Domain *D,LoadList *LL,int s,int iteration)
{
   int l,ii,jj,i,j,istart,iend,jstart,jend,flag,cnt,n,nn,m,numMode;
   int minZSub,minYSub,indexJ;
   double dt,dz,dr,dZ,dR,lambda,weight,charge,weightCoef,delT,coefEr,unitE,invGam2;
	double x,y,r,rr,z,phi,vz,px,py,pz,nenergy,positionZ,positionR,testR,posZ;
	double gamma,gamma0,dGam,delGam,sigGam,emitR,gammaR,distanceR,ne,ne0,rPrime,beta0;
	double sigR,sigRPrime,xPrime,yPrime,totalCnt,phase,density,sinPhi,cosPhi,dPhi;
	double numberRZ,numPhi,numberInCell;
   int myrank,nTasks,rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);	
   Particle **particle;
   particle=D->particle;

   ptclList *New,*p;   

   numMode=D->numMode;
   istart=D->istart; iend=D->iend;
   jstart=D->jstart; jend=D->jend;
   minZSub=D->minXSub;   minYSub=D->minYSub;
   numPhi=LL->numberPhi; numberRZ=LL->numberRZ;
   dZ=D->dz*D->lambda;   dR=D->dr*D->lambda;
	dz=D->dz; 	          dr=D->dr;		dt=D->dt;
   lambda=D->lambda;
		
   charge=LL->charge;

   gamma0=LL->energy/mc2;
   dGam=LL->spread*gamma0;
   emitR=LL->emitR/gamma0;
   gammaR=(1+LL->alphaR*LL->alphaR)/LL->betaR;

   sigR=sqrt(emitR/gammaR)/D->lambda;
   sigRPrime=sqrt(emitR*gammaR);

   distanceR=LL->alphaR/gammaR/D->lambda+LL->focalL/D->lambda;
   vz=sqrt(gamma0*gamma0-1.0)/gamma0;	
   if(vz==0.0)  delT=0.0; 
   else         delT=distanceR/vz;                  //normalized

   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;	

	density=LL->density;

   srand((iteration+1)*myrank);
   numberInCell=numPhi*numberRZ;	
   weight=1.0/numberInCell;

   //position define
   double sum,coef,v[4],tmpMin[nTasks],tmpMax[nTasks],lowGam,upGam,intervalGam;
   double minGam,maxGam;
   const gsl_rng_type * T;

   gsl_rng *ran;
   gsl_rng_env_setup();
   T = gsl_rng_default;
   ran = gsl_rng_alloc(T);
   gsl_qrng *q1=gsl_qrng_alloc(gsl_qrng_reversehalton,4);
   lowGam=1e10; upGam=-1e10;

   //position define     
   for(i=istart; i<iend; i++) 
	{
		z=(i+D->minXSub-istart)*dZ;
	 	ne=0.0;
    	if(LL->gaussMode==OFF) {
      	for(l=0; l<LL->xnodes-1; l++) {
        		posZ=(double)(i+D->minXSub-istart+D->minXDomain);
        		if(posZ>=LL->xpoint[l] && posZ<LL->xpoint[l+1]) {
            	ne=(LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posZ-LL->xpoint[l])+LL->xn[l];
         	} else ;
       	}
	  	} else if(LL->gaussMode==ON) {
	   	posZ=LL->posZ;
	    	phase=pow((z-posZ)/LL->sigZ,LL->gaussPower);
			ne=exp(-phase*0.5);
	  	} else ;

     //energy chirp
     	nenergy=1.0+(z-LL->posZ)*LL->eChirp/LL->energy;

		totalCnt=ne*sigR*sigR/dr/dr*numberInCell*2.0*M_PI;
		cnt=0;
		while(cnt<totalCnt)  {
			 flag=0;
			 while(flag==0) {
				  gsl_qrng_get(q1,v);
				  sum=0.0;
				  for(ii=0; ii<4; ii++) sum+=v[ii];
				    if(sum==0.0) flag=0; else flag=1 ;
		    }

			 coef=sqrt(-2.0*log(v[0]));
			 x=coef*cos(2*M_PI*v[1]);
			 xPrime=coef*sin(2*M_PI*v[1]);
			 x*=sigR;
			 xPrime*=sigRPrime;
			 
			 coef=sqrt(-2.0*log(v[2]));
			 y=coef*cos(2*M_PI*v[3]);
			 yPrime=coef*sin(2*M_PI*v[3]);
			 y*=sigR;
			 yPrime*=sigRPrime;

          sigGam=dGam*0.5*1.201122;    // 1.201122=1/sqrt(ln(2))
          delGam=gaussian_dist(sigGam);
          gamma=gamma0*nenergy+delGam;
          if(gamma<=lowGam) lowGam=gamma; else ;
	       if(gamma>upGam) upGam=gamma; else ;				 

			 pz=sqrt((gamma*gamma-1.0)/(1.0+2*rPrime*rPrime));
   	    px=xPrime*pz;
          py=yPrime*pz;

          x-=delT*px/gamma;
		    y-=delT*py/gamma;
			 r=sqrt(x*x+y*y)/dr;
					 		
			 indexJ=((int)r)+jstart;
			 if(indexJ>=jstart && indexJ<jend) {
				  New = (ptclList *)malloc(sizeof(ptclList)); 
              New->next = particle[i][indexJ].head[s]->pt;
              particle[i][indexJ].head[s]->pt = New;
 	
		        positionZ=randomValue(1.0);
              New->z = positionZ;
   	        New->oldZ= i+positionZ-pz/gamma*dt;
         	  New->x=x/dr;
              New->y=y/dr;
              New->oldX=x/dr-px/gamma*dt;
              New->oldY=y/dr-py/gamma*dt;
	           New->weight=weight;
   	        New->charge=charge;	
      	     New->Ez=New->Ex=New->Ey=0.0;
         	  New->Bz=New->Bx=New->By=0.0;            	     
              New->pz=pz;      New->px=px;       New->py=py;
              D->index+=1;
	           New->index=D->index;            
   	        New->core=myrank;            
		     } else ;
			  cnt++;
		}
 	}			//End of for(i)
   gsl_qrng_free(q1);
   gsl_rng_free(ran);

   if(fabs(lowGam)==1e10) lowGam=0.0; else ;
   if(fabs(upGam)==1e10) upGam=0.0; else ;
   MPI_Gather(&lowGam,1,MPI_DOUBLE,tmpMin,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(tmpMin,nTasks,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Gather(&upGam,1,MPI_DOUBLE,tmpMax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(tmpMax,nTasks,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);
							
   lowGam=tmpMin[0]; upGam=tmpMax[0];
   for(rank=1; rank<nTasks; rank++) {
     if(lowGam>=tmpMin[rank]) lowGam=tmpMin[rank]; else ;
     if(upGam<tmpMax[rank])   upGam=tmpMax[rank]; else ;
   }
   intervalGam=(upGam-lowGam)/(LL->numGam*1.0);
															 
   if(myrank==0) {
     printf("minGamma=%g, maxGamma=%g\n",lowGam,upGam);
   } else ;	

   double **den,**Er,**Ez,**Bp,Wr[2],Wz[2];
	double rho0,alpha,invR,factor[2],tmp;
	int index;

   den = (double **)malloc((D->nxSub+5)*sizeof(double *));
   Er = (double **)malloc((D->nxSub+5)*sizeof(double *));
   Ez = (double **)malloc((D->nxSub+5)*sizeof(double *));
   Bp = (double **)malloc((D->nxSub+5)*sizeof(double *));
   for(i=0; i<D->nxSub+5; i++)   {
     den[i] = (double *)malloc((D->nySub+5)*sizeof(double ));
     Er[i] = (double *)malloc((D->nySub+5)*sizeof(double ));
     Ez[i] = (double *)malloc((D->nySub+5)*sizeof(double ));
     Bp[i] = (double *)malloc((D->nySub+5)*sizeof(double ));
	}
   for(i=0; i<D->nxSub+5; i++) 
     for(j=0; j<D->nySub+5; j++)   {
		  den[i][j]=0.0;
		  Er[i][j]=0.0;
		  Ez[i][j]=0.0;
		  Bp[i][j]=0.0;
	  }

//   rho0=LL->density/LL->criticalDensity;
   rho0=LL->density;
   coefEr=-eCharge/eps0;	
   unitE=eCharge/eMass/D->omega/velocityC;
   alpha=2.0;
   for(n=0; n<LL->numGam; n++) {
     minGam=lowGam+n*intervalGam;
     maxGam=minGam+intervalGam;
     gamma0=(minGam+maxGam)*0.5;

     for(i=0; i<iend+3; i++) 
       for(j=0; j<jend+3; j++)   {
		    den[i][j]=0.0;
          Er[i][j]=0.0;
		    Ez[i][j]=0.0;
		    Bp[i][j]=0.0;
	    }

     for(i=istart; i<iend; i++)
       for(j=jstart+1; j<jend; j++) {
         p=particle[i][j].head[s]->pt;
         while(p) {
           weight=p->weight*rho0*p->charge;
           z=p->z;   x=p->x;   y=p->y;
           r=sqrt(x*x+y*y);   invR=1.0/r;
           index=j-jstart;
           //Wr[0]=((index+1)*(index+1)-r*r)/(2.0*index+1.0);
           Wr[0]=((index+r)*(1-(r-index)))/(2.0*r);  Wr[1]=1.0-Wr[0];
           Wz[1]=z-(int)(z);                         Wz[0]=1.0-Wz[1];
           factor[0]=weight/(2.0*(index));
           factor[1]=weight/(2.0*(index+1.0));

           pz=p->pz; px=p->px; py=p->py;
           gamma=sqrt(1.0+px*px+py*py+pz*pz);
           if(gamma>=minGam && gamma<maxGam) {
	          for(ii=0; ii<2; ii++)
	            for(jj=0; jj<2; jj++) {
                 tmp=Wr[jj]*Wz[ii]*factor[jj];
			        den[i+ii][j+jj]+=tmp;
	            }			  
			  } else ;
 
           p=p->next;
         }
	    }     //End of for(i,j)

     for(i=istart; i<iend; i++)
       for(j=jstart; j<jstart+1; j++) {
         p=particle[i][j].head[s]->pt;
         while(p) {
           weight=p->weight*rho0*p->charge;
           z=p->z;   x=p->x;   y=p->y;
           r=sqrt(x*x+y*y);   invR=1.0/r;
           index=j-jstart;
           if(r<0.5) Wr[0]=(r*r-r+0.5)/(2.0*r*r+0.5);
           else      Wr[0]=(1.0-r)*0.5;
           Wr[1]=1.0-Wr[0];
           Wz[1]=z-(int)(z);              Wz[0]=1.0-Wz[1];

			  pz=p->pz; px=p->px; py=p->py;
           gamma=sqrt(1.0+px*px+py*py+pz*pz);
           if(gamma>=minGam && gamma<maxGam) {
             factor[0]=weight/(0.25);
             factor[1]=weight/(2.0);
	          for(ii=0; ii<2; ii++)
	            for(jj=0; jj<2; jj++) {
                 tmp=Wr[jj]*Wz[ii]*factor[jj];
			        den[i+ii][j+jj]+=tmp;
	            }			  
			  } else ;
 
           p=p->next;
         }
	    }     //End of for(i,j)

     if(nTasks>1)  {     
       MPI_TransferBeamDen_Xplus(D,den,D->nySub+5,3);
       MPI_TransferBeamDen_Xminus(D,den,D->nySub+5,3);
	  } else ;
     if(nTasks>1)  {     
       MPI_TransferBeamField_Xplus(D,den,D->nySub+5,3);
       MPI_TransferBeamField_Xminus(D,den,D->nySub+5,3);
	  } else ;

     // Cal Er
     for(i=1; i<iend+1; i++) {
		 j=jstart;
		   Er[i][j]=0.0;  
       for(j=jstart+1; j<jend+3; j++) {
		   r=(j-jstart)*dR;
		   sum=0.0;
         for(jj=jstart+1; jj<j; jj++) {
		     rr=(jj-jstart)*dR;
			  sum+=den[i][jj]*rr*dR;
			}
			Er[i][j]=sum*coefEr/r;
		 }
	  }
     if(nTasks>1)  {     
       MPI_TransferBeamField_Xplus(D,Er,D->nySub+5,3);
       MPI_TransferBeamField_Xminus(D,Er,D->nySub+5,3);
	  } else ;

     // Cal Ez
     invGam2=1.0/gamma0/gamma0;
     for(i=istart; i<iend; i++) 
       for(j=jstart; j<jend+2; j++) {
         sum=0.0;
         for(jj=j; jj<jend; jj++)
           sum+=Er[i+1][jj]-Er[i-1][jj];
         Ez[i][j]=sum*invGam2/(2.0*dZ)*dR;
		 }
     if(nTasks>1)  {     
       MPI_TransferBeamField_Xplus(D,Ez,D->nySub+5,3);
       MPI_TransferBeamField_Xminus(D,Ez,D->nySub+5,3);
	  } else ;

     // Cal Bphi
	  beta0=sqrt(1.0-invGam2);
     for(i=0; i<iend+3; i++) 
       for(j=0; j<jend+3; j++)
         Bp[i][j]=Er[i][j]*beta0;

	  int shift=-0;
     switch(D->fieldType)  {
     case Yee :
     case NoCherenkov :
       for(i=istart-1; i<iend+3; i++) 
         for(j=jstart; j<jend+3; j++) {
           D->EzR[0][i][j]+=Ez[i][j];
           D->ErR[0][i][j]-=Er[i][j];
           D->BpR[0][i][j]-=Bp[i][j];
//           D->EzR[0][i][j]+=0.5*(Ez[i-1][j]+Ez[i][j])*unitE;
//           D->ErR[0][i][j]-=0.5*(Er[i-1][j]+Er[i-1][j+1])*unitE;
//           D->BpR[0][i][j]-=0.5*(Bp[i-1][j]+Bp[i-1][j+1])*unitE;
	 }
       break;
     case Split :
       
       for(i=istart-1; i<iend+3; i++) 
         for(j=jstart; j<jend+3; j++) {
           D->EzR[0][i][j]+=Ez[i][j]*unitE;	
           D->EzNowR[0][i][j]+=Ez[i][j]*unitE;	
           D->PrR[0][i][j]+=(Er[i][j]+Bp[i][j])*unitE;
           D->PlR[0][i][j]+=(Er[i][j]-Bp[i][j])*unitE;
//           D->PrR[0][i][j]-=0.5*(Er[i][j]+Bp[i][j]+Er[i][j+1]+Bp[i][j+1])*unitE;
//           D->PlR[0][i][j]-=0.5*(Er[i][j]-Bp[i][j]+Er[i][j+1]-Bp[i][j+1])*unitE;
	 		}
       break;
     }	


   }

   LL=D->loadList; s=0;
   while(LL->next)      {
      if(LL->type==Beam)
         solveCharge(D,LL,D->RhoNoPairR,D->RhoNoPairI,istart,iend,jstart,jend,s,1.0);
      else ;
      LL=LL->next; s++;
   }
   D->shareF[0]=D->RhoNoPairR;
   MPI_TransferJNew_Xplus(D,1,D->nySub+5,3);
   MPI_TransferJNew_Xminus(D,1,D->nySub+5,3);
   if(D->Period==ON) { MPI_TransferJNew_Period_X(D,1,D->nySub+5,3); } else ;
   D->RhoNoPairR=D->shareF[0];

   D->shareF[0]=D->RhoNoPairR;
   MPI_TransferFNew_Xplus(D,1,D->nySub+5,3);
   MPI_TransferFNew_Xminus(D,1,D->nySub+5,3);
   if(D->Period==ON) { MPI_TransferFNew_Period_X(D,1,D->nySub+5,3); } else ;
   D->RhoNoPairR=D->shareF[0];

   switch(D->fieldType)  {
   case Yee :
   case NoCherenkov :
     for(i=istart-1; i<iend+3; i++) 
       for(j=jstart; j<jend+3; j++) {
         D->EzR[0][i][j]+=Ez[i][j];
         D->ErR[0][i][j]-=Er[i][j];
         D->BpR[0][i][j]-=Bp[i][j];
	    }
       break;
   case Split :       
     m=0;
       for(i=istart-1; i<iend+2; i++) 
         for(j=jstart+1; j<jend+3; j++) {
           r=j-jstart;
           D->PrR[m][i][j]=2*r/(r+0.5)*(
                             M_PI*dr*(D->RhoNoPairR[m][i+1][j]+D->RhoNoPairR[m][i][j])
                             -dr/dz*(D->EzR[m][i+1][j]-D->EzR[m][i][j])
                           ) +(r-0.5)/(r+0.5)*D->PrR[m][i][j-1];
	 		}
       j=jstart;
       for(i=istart-1; i<iend+2; i++) { 
           D->PrR[m][i][j]=0.5*D->PrR[m][i][j];
	    }
     break;
   }	
   D->shareF[0]=D->PrR;
   MPI_TransferFNew_Xplus(D,1,D->nySub+5,3);
   MPI_TransferFNew_Xminus(D,1,D->nySub+5,3);
   if(D->Period==ON) { MPI_TransferFNew_Period_X(D,1,D->nySub+5,3); } else ;
   D->PrR=D->shareF[0];
   
//lala

   //FILE *out;
	//out=fopen("density","w");
   //for(i=istart; i<iend; i++) {
   //  for(j=jstart; j<jend; j++) {
	//	  z=(i-istart+D->minXDomain)*dZ;
	//	  r=(j-jstart+D->minYDomain)*dR;
	//	  fprintf(out,"%g %g %g %g %g\n",z,r,den[i][j],Er[i][j]*unitE,Ez[i][j]*unitE);
	//  }
	//  fprintf(out,"\n");
   //}
	//fclose(out);



   for(i=0; i<D->nxSub+5; i++) {
     free(den[i]); 
     free(Er[i]); 
     free(Ez[i]); 
     free(Bp[i]); 
	}
   free(den); 
   free(Er); 
   free(Ez); 
   free(Bp); 
   

}

double gaussian_dist(double sig)
{
  double r,prob,v,random;
  int intRand,randRange=1e5;

  r=1.0;
  prob=0.0;
  while (r>prob)  {
    intRand = rand() % randRange;
    r = ((double)intRand)/randRange;
    intRand = rand() % randRange;
    random = ((double)intRand)/randRange;
    v = 4.0*(random-0.5);
    prob=exp(-v*v);
  }
  return sig*v;
}

void MPI_TransferBeamField_Xplus(Domain *D,double **f1,int ny,int share)
{
  int i,j,num,start;
  int istart,iend,jstart,jend;  
  int myrank, nTasks;
  double *data;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart; iend=D->iend;
  jstart=D->istart; jend=D->iend;

  num=ny*(share-1);
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores
  start=0;
  for(i=1; i<share; i++) {
    for(j=0; j<ny; j++) data[j+start]=f1[iend-i][j]; start+=ny;
  }

  if(myrank%2==1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);	  
    start=0;
    for(i=1; i<share; i++) {
      for(j=0; j<ny; j++)  f1[istart-i][j]=data[j+start];  start+=ny;
    }
  }
  else if(myrank%2==0 && myrank!=nTasks-1)
    MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);	  

  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores
  start=0;
  for(i=1; i<share; i++) {
    for(j=0; j<ny; j++) data[j+start]=f1[iend-i][j]; start+=ny;
  }

  if(myrank%2==0 && myrank!=0)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);	  
    start=0;
    for(i=1; i<share; i++) {
      for(j=0; j<ny; j++)  f1[istart-i][j]=data[j+start];  start+=ny;
    }
  }
  else if(myrank%2==1 && myrank!=nTasks-1)
    MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);	  
  MPI_Barrier(MPI_COMM_WORLD);

  free(data);
}

void MPI_TransferBeamField_Xminus(Domain *D,double **f1,int ny,int share)
{
  int i,j,num,start;
  int istart,iend,jstart,jend;  
  int myrank, nTasks;
  double *data;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart; iend=D->iend;
  jstart=D->istart; jend=D->iend;

  num=ny*share;
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores
  start=0;
  for(i=0; i<share; i++) {
    for(j=0; j<ny; j++) data[j+start]=f1[istart+i][j]; start+=ny;
  }

  if(myrank%2==0 && myrank!=nTasks-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);	  
    start=0;
    for(i=0; i<share; i++) {
      for(j=0; j<ny; j++)  f1[iend+i][j]=data[j+start];  start+=ny;
    }
  }
  else if(myrank%2==1)
    MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);	  

  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores
  start=0;
  for(i=0; i<share; i++) {
    for(j=0; j<ny; j++) data[j+start]=f1[istart+i][j]; start+=ny;
  }

  if(myrank%2==1 && myrank!=nTasks-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);	  
    start=0;
    for(i=0; i<share; i++) {
      for(j=0; j<ny; j++)  f1[iend+i][j]=data[j+start];  start+=ny;
    }
  }
  else if(myrank%2==0 && myrank!=0)
    MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);	  
  MPI_Barrier(MPI_COMM_WORLD);

  free(data);
}


void MPI_TransferBeamDen_Xplus(Domain *D,double **f1,int ny,int share)
{
  int i,j,num,start;
  int istart,iend,jstart,jend;  
  int myrank, nTasks;
  double *data;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart; iend=D->iend;
  jstart=D->istart; jend=D->iend;

  num=ny*share;
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores
  start=0;
  for(i=0; i<share; i++) {
    for(j=0; j<ny; j++) data[j+start]=f1[iend+i][j]; start+=ny;
  }

  if(myrank%2==1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);	  
    start=0;
    for(i=0; i<share; i++) {
      for(j=0; j<ny; j++)  f1[istart+i][j]+=data[j+start];  start+=ny;
    }
  }
  else if(myrank%2==0 && myrank!=nTasks-1)
    MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);	  

  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores
  start=0;
  for(i=0; i<share; i++) {
    for(j=0; j<ny; j++) data[j+start]=f1[iend+i][j]; start+=ny;
  }

  if(myrank%2==0 && myrank!=0)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);	  
    start=0;
    for(i=0; i<share; i++) {
      for(j=0; j<ny; j++)  f1[istart+i][j]+=data[j+start];  start+=ny;
    }
  }
  else if(myrank%2==1 && myrank!=nTasks-1)
    MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);	  
  MPI_Barrier(MPI_COMM_WORLD);

  free(data);
}

void MPI_TransferBeamDen_Xminus(Domain *D,double **f1,int ny,int share)
{
  int i,j,num,start;
  int istart,iend,jstart,jend;  
  int myrank, nTasks;
  double *data;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart; iend=D->iend;
  jstart=D->istart; jend=D->iend;

  num=ny*(share-1);
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores
  start=0;
  for(i=1; i<share; i++) {
    for(j=0; j<ny; j++) data[j+start]=f1[istart-i][j]; start+=ny;
  }

  if(myrank%2==0 && myrank!=nTasks-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);	  
    start=0;
    for(i=1; i<share; i++) {
      for(j=0; j<ny; j++)  f1[iend-i][j]+=data[j+start];  start+=ny;
    }
  }
  else if(myrank%2==1)
    MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);	  

  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores
  start=0;
  for(i=1; i<share; i++) {
    for(j=0; j<ny; j++) data[j+start]=f1[istart-i][j]; start+=ny;
  }

  if(myrank%2==1 && myrank!=nTasks-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);	  
    start=0;
    for(i=1; i<share; i++) {
      for(j=0; j<ny; j++)  f1[iend-i][j]+=data[j+start];  start+=ny;
    }
  }
  else if(myrank%2==0 && myrank!=0)
    MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);	  
  MPI_Barrier(MPI_COMM_WORLD);

  free(data);
}

