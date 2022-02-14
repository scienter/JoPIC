#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_qrng.h>

double randomValue(double beta);
double maxwellianVelocity(double temperature,double mass);
void loadMovingPlasma_crystal(Domain *D,LoadList *LL,int s);
double applyFunctionX(int mode,double centerX,double x,double gaussCoefX,double polyCoefX);
double applyFunctionYZ(int mode,double centerY,double y,double centerZ,double z,double gaussCoefYZ,double polyCoefYZ);
void random1D_sobol(double *x,gsl_qrng *q);
void random2D_sobol(double *x,double *y,gsl_qrng *q);
void random3D_sobol(double *x,double *y,double *z,gsl_qrng *q);



void assignDefParticle(Domain *D)
{
  int n,cnt,flag;
  double tmp,tmpX,tmpY,tmpZ,tmpR,factY,factZ,randX,randY,randZ;
  LoadList *LL,*prevL;
  DefPtcl *p,*prevP;

  factY=factZ=0.0;
  if(D->dimension>1) factY=1.0; else;
  if(D->dimension>2) factZ=1.0; else;

  LL=D->loadList;
  while(LL->next)      {
    if(LL->type==Defined) {
      if(LL->pair==OFF) {
        p=LL->def;
        while(p) {
          if(p->flag==OFF) {
            if(p->xPos<=LL->createGuard) {
              cnt=p->numDefPtcls=LL->numDefPtcls;           
              p->define=(double *)malloc((cnt*3)*sizeof(double ));
              gsl_qrng *q3 = gsl_qrng_alloc (gsl_qrng_sobol,3);
              for(n=0; n<cnt; n++)  {
                flag=0;
                while(flag==0) {
//                  random2D_sobol(&randX,&randY,q3); randZ=0.0;
                  random3D_sobol(&randX,&randY,&randZ,q3);
                  tmpX=-LL->RDef+randX*2.0*LL->RDef;
                  tmpY=(-LL->RDef+randY*2.0*LL->RDef)*factY;
                  tmpZ=(-LL->RDef+randZ*2.0*LL->RDef)*factZ;
//                  tmp=(double)(randomValue(1.0));
//                  tmpX=-LL->RDef+tmp*2.0*LL->RDef;
//                  tmp=(double)(randomValue(1.0));
//                  tmpY=(-LL->RDef+tmp*2.0*LL->RDef)*factY;
//                  tmp=(double)(randomValue(1.0));
//                  tmpZ=(-LL->RDef+tmp*2.0*LL->RDef)*factZ;
                  tmpR=sqrt(tmpX*tmpX+tmpY*tmpY+tmpZ*tmpZ);
                  if(tmpR<=LL->RDef) {
                    flag=1;
                    tmpX=p->xPos+tmpX;
                    tmpY=p->yPos+tmpY;
                    tmpZ=p->zPos+tmpZ;
                  } else  flag=0;
                }
                p->define[n]=tmpX;
                p->define[n+cnt]=tmpY;
                p->define[n+2*cnt]=tmpZ;
              }	//End of for(n<cnt)
              p->flag=ON;  
            } else ;
          } else ;	//End of if(flag==Off)
          p=p->next;
        }
      } else { 	//if pair==ON
        p=LL->def; prevP=prevL->def;
        while(p) {
          if(p->flag==OFF) {
            if(p->xPos<=LL->createGuard) {
              cnt=p->numDefPtcls=prevL->numDefPtcls;           
              p->define=(double *)malloc((cnt*3)*sizeof(double ));
              for(n=0; n<cnt; n++)  {
                p->define[n]=prevP->define[n];
                p->define[n+cnt]=prevP->define[n+cnt];
                p->define[n+2*cnt]=prevP->define[n+2*cnt];
              }
              p->flag=ON;  
            } else ;
          } else ;	//End of if(flag==Off)
          p=p->next;
          prevP=prevP->next;	
        } 
      }
    } else ; 		//End of if(LL->type==Defined) 
    prevL=LL;
    LL=LL->next;
  }
}

void deleteDefParticle(Domain *D)
{
  int cnt;
  LoadList *LL;
  DefPtcl *p,*prev;

  LL=D->loadList;
  while(LL->next)      {
    if(LL->type==Defined) {
      p=LL->def;
      cnt=1;
      while(p) {
        if(cnt==1) prev=p; else ;
  
        if(p->xPos<=LL->deleteGuard) {
          if(cnt==1)  {
            LL->def=p->next;
            p->next=NULL;
            free(p->define);
            free(p);
            p=LL->def;
            cnt=1;
          } else {
            prev->next=p->next;
            p->next=NULL;
            free(p->define);
            free(p);
            p=prev->next;
          } 
        }
        else  {
          prev=p;
          p=p->next;
          cnt++;
        }
      }
    }	else ;	//End of if(LL->type==Defined)
    LL=LL->next;
  }
}

void saveDefParticle(Domain *D,int iteration)
{
  int n,cnt,num;
  double x,y,z;
  char name[100];
  FILE *out;

  LoadList *LL;
  DefPtcl *p;
  LL=D->loadList;
  num=0;
  while(LL->next)      {
    sprintf(name,"def%d_%d",iteration,num);
    out = fopen(name,"w");    

    p=LL->def;
    while(p) {
      cnt=p->numDefPtcls;
      for(n=0; n<cnt; n++)  {
        x=p->define[n]*D->lambda;
        y=p->define[n+cnt]*D->lambda;
        z=p->define[n+2*cnt]*D->lambda;
        fprintf(out,"%g %g %g\n",x,y,z);
      }	//End of for(n<cnt)
      p=p->next;
    }
    fclose(out);
    num++;
    LL=LL->next;
  }
}
