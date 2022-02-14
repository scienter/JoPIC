#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "mpi.h"

void interpolation_Split_1st(Domain *D,External *Ext);
void interpolation_Split_2nd(Domain *D,External *Ext);
void interpolation_Yee_Pukhov_1st(Domain *D,External *Ext);
void interpolation_Yee_Pukhov_2nd(Domain *D,External *Ext);
double secondInterpolation(double ***F,double ***Fxyz,double **Fxy,double *Fx,int i,int j,int k,double x,double y,double z,
                          int jjN,int kkN,double coef2D,double coef3D,int grid2D,int grid3D);

void interpolation(Domain D,External *Ext)
{
//  switch((D.fieldType-1)*6+(D.interpolationType-1)*3+D.dimension) {  
  switch( (D.fieldType-1)*2+D.interpolationType) {  

  case ( (Split-1)*2+FIRST ) :
    interpolation_Split_1st(&D,Ext);
    break;
  case ( (Split-1)*2+SECOND ) :
    interpolation_Split_2nd(&D,Ext);
    break;

  //Yee, Pukhov
  case ( (Yee-1)*2+FIRST ) :
  case ( (Pukhov-1)*2+FIRST ) :
    interpolation_Yee_Pukhov_1st(&D,Ext);
    break;
  case ( (Yee-1)*2+SECOND ) :
  case ( (Pukhov-1)*2+SECOND ) :
    interpolation_Yee_Pukhov_2nd(&D,Ext);
    break;

  default :
    printf("In interpolation, what interpolationType(%d)? and what dimension(%d)?\n",D.interpolationType,D.dimension);
  }
}


void interpolation_Split_1st(Domain *D,External *Ext)
{
  int i,j,k,i1,j1,k1,istart,iend,jstart,jend,kstart,kend,s,cnt,ii,jj,kk;
  double E1,Pr,Pl,B1,Sr,Sl,extE1,extE2,extE3,extB1,extB2,extB3,x,y,z,x1,y1,z1;
	double wx1[2],wy1[2],wz1[2],wx2[2],wy2[2],wz2[2];
  int iiN,jjN,kkN,grid2D,grid3D,effect2D,effect3D;
  ptclList *p;
  int myrank, nprocs;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  kstart=D->kstart;   kend=D->kend;
  Particle ***particle;
  particle=D->particle;
  extE1=Ext->E1;
  extE2=Ext->E2;
  extE3=Ext->E3;
  extB1=Ext->B1;
  extB2=Ext->B2;
  extB3=Ext->B3;

  iiN=2; jjN=1; kkN=1;
  effect2D=effect3D=0.0;
  grid2D=1; grid3D=1;
  if(D->dimension>1) { jjN=2; effect2D=1; grid2D=0; } else ;
  if(D->dimension>2) { kkN=2; effect3D=1; grid3D=0; } else ; 

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      for(k=kstart; k<kend; k++)
      {
        for(s=0; s<D->nSpecies; s++)
        {
          p=particle[i][j][k].head[s]->pt;
          while(p)  {
            x=p->x;  y=p->y;  z=p->z;
            i1=(int)(i+x+0.5);       j1=(int)(j+y+0.5);       k1=(int)(k+z+0.5);
            x1=x+0.5-((int)(x+0.5)); y1=y+0.5-((int)(y+0.5)); z1=z+0.5-((int)(z+0.5));
            wx1[1]=x;		          wx1[0]=1.0-wx1[1];
            wy1[1]=y*effect2D;		wy1[0]=1.0-wy1[1];
            wz1[1]=z*effect3D;		wz1[0]=1.0-wz1[1];
            wx2[1]=x1;		        wx2[0]=1.0-wx2[1];
            wy2[1]=y1*effect2D;		wy2[0]=1.0-wy2[1];
            wz2[1]=z1*effect3D;		wz2[0]=1.0-wz2[1];

			      B1=E1=Pr=Pl=Sr=Sl=0.0;
            for(ii=0; ii<iiN; ii++)
              for(jj=0; jj<jjN; jj++)
                for(kk=0; kk<kkN; kk++) {
                  E1+=wx2[ii]*wy1[jj]*wz1[kk]*D->ExNow[i1-1+ii][j+jj][k+kk];
		              Pr+=wx2[ii]*wy2[jj]*wz1[kk]*D->Pr[i1-1+ii][j1-1+jj+grid2D][k+kk];
		              Pl+=wx2[ii]*wy2[jj]*wz1[kk]*D->Pl[i1-1+ii][j1-1+jj+grid2D][k+kk];
		              B1+=wx2[ii]*wy2[jj]*wz2[kk]*D->BxNow[i1-1+ii][j1-1+jj+grid2D][k1-1+kk+grid3D];
		              Sr+=wx2[ii]*wy1[jj]*wz2[kk]*D->Sr[i1-1+ii][j+jj][k1-1+kk+grid3D];
		              Sl+=wx2[ii]*wy1[jj]*wz2[kk]*D->Sl[i1-1+ii][j+jj][k1-1+kk+grid3D];
				        }
            p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
            p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;
             p=p->next;
          }
        }		//for(s)        

        // tracking particle
        s=0;
        p=D->track[i][j][k].head[s]->pt;
        while(p)  {
          x=p->x;  y=p->y;  z=p->z;
          i1=(int)(i+x+0.5);       j1=(int)(j+y+0.5);       k1=(int)(k+z+0.5);
          x1=x+0.5-((int)(x+0.5)); y1=y+0.5-((int)(y+0.5)); z1=z+0.5-((int)(z+0.5));
          wx1[1]=x;		wx1[0]=1.0-wx1[1];
          wy1[1]=y*effect2D;		wy1[0]=1.0-wy1[1];
          wz1[1]=z*effect3D;		wz1[0]=1.0-wz1[1];
          wx2[1]=x1;	          wx2[0]=1.0-wx2[1];
          wy2[1]=y1*effect2D;	  wy2[0]=1.0-wy2[1];
          wz2[1]=z1*effect3D;	  wz2[0]=1.0-wz2[1];
          for(ii=0; ii<iiN; ii++)
            for(jj=0; jj<jjN; jj++)
              for(kk=0; kk<kkN; kk++) {
		            B1+=wx2[ii]*wy2[jj]*wz2[kk]*D->BxNow[i1-1+ii][j1-1+jj+grid2D][k1-1+kk+grid3D];
                E1+=wx2[ii]*wy1[jj]*wz1[kk]*D->ExNow[i1-1+ii][j+jj][k+kk];
		            Pr+=wx2[ii]*wy2[jj]*wz1[kk]*D->Pr[i1-1+ii][j1-1+jj+grid2D][k+kk];
		            Pl+=wx2[ii]*wy2[jj]*wz1[kk]*D->Pl[i1-1+ii][j1-1+jj+grid2D][k+kk];
		            Sr+=wx2[ii]*wy1[jj]*wz2[kk]*D->Sr[i1-1+ii][j+jj][k1-1+kk+grid3D];
		            Sl+=wx2[ii]*wy1[jj]*wz2[kk]*D->Sl[i1-1+ii][j+jj][k1-1+kk+grid3D];
			       }
          p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
          p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;
          p=p->next;
        }
        
			}		   //for(i,j)

}

void interpolation_Yee_Pukhov_1st(Domain *D,External *Ext)
{
  int i,j,k,i1,j1,k1,ii,jj,kk,istart,iend,jstart,jend,kstart,kend,s;
  double Bx,By,Bz,Ex,Ey,Ez,extE1,extE2,extE3,extB1,extB2,extB3;
  double x,y,z,x1,y1,z1;
  double wx1[2],wy1[2],wz1[2],wx2[2],wy2[2],wz2[2];
  ptclList *p;
  int myrank, nprocs;  
  int iiN,jjN,kkN,grid2D,grid3D,effect2D,effect3D;  

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  kstart=D->kstart;   kend=D->kend;

  Particle ***particle;
  particle=D->particle;

  extE1=Ext->E1;   extE2=Ext->E2;
  extE3=Ext->E3;   extB1=Ext->B1;
  extB2=Ext->B2;   extB3=Ext->B3;

  iiN=2; jjN=1; kkN=1;
  effect2D=effect3D=0.0;
  grid2D=1; grid3D=1;
  if(D->dimension>1) { jjN=2; effect2D=1; grid2D=0; } else ;
  if(D->dimension>2) { kkN=2; effect3D=1; grid3D=0; } else ; 

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      for(k=kstart; k<kend; k++)
      {
        for(s=0; s<D->nSpecies; s++)
        {
          p=particle[i][j][k].head[s]->pt;
          while(p)  {
            x=p->x;  y=p->y;  z=p->z;
            i1=(int)(i+x+0.5);       j1=(int)(j+y+0.5);       k1=(int)(k+z+0.5);
            x1=x+0.5-((int)(x+0.5)); y1=y+0.5-((int)(y+0.5)); z1=z+0.5-((int)(z+0.5));
            wx1[1]=x;		          wx1[0]=1.0-wx1[1];
            wy1[1]=y*effect2D;		wy1[0]=1.0-wy1[1];
            wz1[1]=z*effect3D;		wz1[0]=1.0-wz1[1];
            wx2[1]=x1;		        wx2[0]=1.0-wx2[1];
            wy2[1]=y1*effect2D;		wy2[0]=1.0-wy2[1];
            wz2[1]=z1*effect3D;		wz2[0]=1.0-wz2[1];


            Bx=0.0; By=0.0; Bz=0.0;
            Ex=0.0; Ey=0.0; Ez=0.0;
            for(ii=0; ii<iiN; ii++)
              for(jj=0; jj<jjN; jj++)
                for(kk=0; kk<kkN; kk++) {
                  Bx+=wx1[ii]*wy2[jj]*wz2[kk]*D->BxNow[i+ii][j1-1+jj+grid2D][k1-1+kk+grid3D];
                  By+=wx2[ii]*wy1[jj]*wz2[kk]*D->ByNow[i1-1+ii][j+jj][k1-1+kk+grid3D];
                  Bz+=wx2[ii]*wy2[jj]*wz1[kk]*D->BzNow[i1-1+ii][j1-1+jj+grid2D][k+kk];
                  Ex+=wx2[ii]*wy1[jj]*wz1[kk]*D->Ex[i1-1+ii][j+jj][k+kk];
                  Ey+=wx1[ii]*wy2[jj]*wz1[kk]*D->Ey[i+ii][j1-1+jj+grid2D][k+kk];
                  Ez+=wx1[ii]*wy1[jj]*wz2[kk]*D->Ez[i+ii][j+jj][k1-1+kk+grid3D];
                }
            p->E1=Ex+extE1; p->E2=Ey+extE2; p->E3=Ez+extE3;
            p->B1=Bx+extB1; p->B2=By+extB2; p->B3=Bz+extB3;

            p=p->next;
          } 
        }	//End of for(s)

        //  tracking particle
        s=0;
          p=D->track[i][j][k].head[s]->pt;
          while(p)  {
            x=p->x;  y=p->y;  z=p->z;
            i1=(int)(i+x+0.5);       j1=(int)(j+y+0.5);       k1=(int)(k+z+0.5);
            x1=x+0.5-((int)(x+0.5)); y1=y+0.5-((int)(y+0.5)); z1=z+0.5-((int)(z+0.5));
            wx1[1]=x;		          wx1[0]=1.0-wx1[1];
            wy1[1]=y*effect2D;		wy1[0]=1.0-wy1[1];
            wz1[1]=z*effect3D;		wz1[0]=1.0-wz1[1];
            wx2[1]=x1;		        wx2[0]=1.0-wx2[1];
            wy2[1]=y1*effect2D;		wy2[0]=1.0-wy2[1];
            wz2[1]=z1*effect3D;		wz2[0]=1.0-wz2[1];


            Bx=0.0; By=0.0; Bz=0.0;
            Ex=0.0; Ey=0.0; Ez=0.0;
            for(ii=0; ii<iiN; ii++)
              for(jj=0; jj<jjN; jj++)
                for(kk=0; kk<kkN; kk++) {
                  Bx+=wx1[ii]*wy2[jj]*wz2[kk]*D->BxNow[i+ii][j1-1+jj+grid2D][k1-1+kk+grid3D];
                  By+=wx2[ii]*wy1[jj]*wz2[kk]*D->ByNow[i1-1+ii][j+jj][k1-1+kk+grid3D];
                  Bz+=wx2[ii]*wy2[jj]*wz1[kk]*D->BzNow[i1-1+ii][j1-1+jj+grid2D][k+kk];
                  Ex+=wx2[ii]*wy1[jj]*wz1[kk]*D->Ex[i1-1+ii][j+jj][k+kk];
                  Ey+=wx1[ii]*wy2[jj]*wz1[kk]*D->Ey[i+ii][j1-1+jj+grid2D][k+kk];
                  Ez+=wx1[ii]*wy1[jj]*wz2[kk]*D->Ez[i+ii][j+jj][k1-1+kk+grid3D];
                }
            p->E1=Ex+extE1; p->E2=Ey+extE2; p->E3=Ez+extE3;
            p->B1=Bx+extB1; p->B2=By+extB2; p->B3=Bz+extB3;

            p=p->next;
          } 
       
      }	//End of for(i,j,k)       
}


void interpolation_Split_2nd(Domain *D,External *Ext)
{
  int i,j,k,i1,j1,k1,istart,iend,jstart,jend,kstart,kend,s;
  double x,y,z,x1,y1,z1,Pr,Pl,Sr,Sl,E1,B1;
  double extE1,extE2,extE3,extB1,extB2,extB3;
  double Fx[3],***Fxyz,**Fxy;
  ptclList *p;
  int myrank, nprocs; 
  int iiN,jjN,kkN,grid2D,grid3D,coef2D,coef3D;     

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  kstart=D->kstart;   kend=D->kend;
  Particle ***particle;
  particle=D->particle;
  extE1=Ext->E1;   extE2=Ext->E2;   extE3=Ext->E3;
  extB1=Ext->B1;   extB2=Ext->B2;   extB3=Ext->B3;

  jjN=1; kkN=1;
  coef2D=coef3D=0.0;
  grid2D=1; grid3D=1;
  if(D->dimension>1) { jjN=3; coef2D=1; grid2D=0; } else ;
  if(D->dimension>2) { kkN=3; coef3D=1; grid3D=0; } else ; 
  Fxyz=(double ***)malloc(3*sizeof(double **));
  Fxy=(double **)malloc(3*sizeof(double *));
  for(i=0; i<3; i++) {
    Fxyz[i]=(double **)malloc(3*sizeof(double *));
    Fxy[i]=(double *)malloc(3*sizeof(double ));
    for(j=0; j<3; j++) 
      Fxyz[i][j]=(double *)malloc(3*sizeof(double ));
  }
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++) { 
      for(k=0; k<3; k++) 
        Fxyz[i][j][k]=0.0;
      Fxy[i][j]=0.0;
    }
    Fx[i]=0.0;
  }

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      for(k=kstart; k<kend; k++)
      {
        for(s=0; s<D->nSpecies; s++)
        {
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {
            x=p->x;  y=p->y;  z=p->z;
            i1=(int)(i+x+0.5);
            j1=(int)(j+y+0.5);
            k1=(int)(k+z+0.5);
            x1=x+0.5;  y1=y+0.5;  z1=z+0.5;

            Pr=secondInterpolation(D->Pr,Fxyz,Fxy,Fx,i,j,k1,x1,y1,z,jjN,kkN,coef2D,coef3D,grid2D,grid3D);
            Pl=secondInterpolation(D->Pl,Fxyz,Fxy,Fx,i,j,k1,x1,y1,z,jjN,kkN,coef2D,coef3D,grid2D,grid3D);
            Sr=secondInterpolation(D->Sr,Fxyz,Fxy,Fx,i,j1,k,x1,y,z1,jjN,kkN,coef2D,coef3D,grid2D,grid3D);
            Sl=secondInterpolation(D->Sl,Fxyz,Fxy,Fx,i,j1,k,x1,y,z1,jjN,kkN,coef2D,coef3D,grid2D,grid3D);
            E1=secondInterpolation(D->ExNow,Fxyz,Fxy,Fx,i,j1,k1,x1,y,z,jjN,kkN,coef2D,coef3D,grid2D,grid3D);
            B1=secondInterpolation(D->BxNow,Fxyz,Fxy,Fx,i,j,k,x1,y1,z1,jjN,kkN,coef2D,coef3D,grid2D,grid3D);

            p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
            p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;
     
            p=p->next;
          }
        }		//for(s)        

        //  tracking particle
        s=0;
        p=D->track[i][j][k].head[s]->pt;
        while(p)    {
          x=p->x;  y=p->y;  z=p->z;
          i1=(int)(i+x+0.5);
          j1=(int)(j+y+0.5);
          k1=(int)(k+z+0.5);
          x1=x+0.5;  y1=y+0.5;  z1=z+0.5;

          Pr=secondInterpolation(D->Pr,Fxyz,Fxy,Fx,i,j,k1,x1,y1,z,jjN,kkN,coef2D,coef3D,grid2D,grid3D);
          Pl=secondInterpolation(D->Pl,Fxyz,Fxy,Fx,i,j,k1,x1,y1,z,jjN,kkN,coef2D,coef3D,grid2D,grid3D);
          Sr=secondInterpolation(D->Sr,Fxyz,Fxy,Fx,i,j1,k,x1,y,z1,jjN,kkN,coef2D,coef3D,grid2D,grid3D);
          Sl=secondInterpolation(D->Sl,Fxyz,Fxy,Fx,i,j1,k,x1,y,z1,jjN,kkN,coef2D,coef3D,grid2D,grid3D);
          E1=secondInterpolation(D->ExNow,Fxyz,Fxy,Fx,i,j1,k1,x1,y,z,jjN,kkN,coef2D,coef3D,grid2D,grid3D);
          B1=secondInterpolation(D->BxNow,Fxyz,Fxy,Fx,i,j,k,x1,y1,z1,jjN,kkN,coef2D,coef3D,grid2D,grid3D);

          p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
          p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;

          p=p->next;
        }

      }		   //for(i,j,k)


  for(i=0; i<3; i++) {
    for(j=0; j<3; j++)
      free(Fxyz[i][j]);
    free(Fxy[i]);
    free(Fxyz[i]);
  }
  free(Fxyz);
  free(Fxy);

}



void interpolation_Yee_Pukhov_2nd(Domain *D,External *Ext)
{
  int i,j,k,ii,jj,kk,i1,j1,k1,istart,iend,jstart,jend,kstart,kend,s;
  double x,y,z,E1,E2,E3,B1,B2,B3;
  double Wx[3],Wy[3],Wz[3],WxC[3],WyC[3],WzC[3];
  double extE1,extE2,extE3,extB1,extB2,extB3;
  ptclList *p;
  int myrank, nprocs;   
  int iiN,jjN,kkN,grid2D,grid3D,effect2D,effect3D;    

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  kstart=D->kstart;   kend=D->kend;

  Particle ***particle;
  particle=D->particle;

  extE1=Ext->E1;   extE2=Ext->E2;   extE3=Ext->E3;
  extB1=Ext->B1;   extB2=Ext->B2;   extB3=Ext->B3;

  iiN=3; jjN=1; kkN=1;
  effect2D=effect3D=0.0;
  grid2D=1; grid3D=1;
  if(D->dimension>1) { jjN=3; effect2D=1; grid2D=0; } else ;
  if(D->dimension>2) { kkN=3; effect3D=1; grid3D=0; } else ;  


  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      for(k=kstart; k<kend; k++)
      {
        for(s=0; s<D->nSpecies; s++)   {
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {
            x=p->x;  y=p->y;  z=p->z;
            WxC[2]=0.5*x*x;
            WxC[1]=0.75-(0.5-x)*(0.5-x);
            WxC[0]=1.0-WxC[1]-WxC[2]; //0.5*(1-x)*(1-x);
            WyC[1]=(0.75-(0.5-y)*(0.5-y))*effect2D;
            WyC[2]=0.5*y*y*effect2D;
            WyC[0]=1.0-WyC[1]-WyC[2]; //0.5*(1-y)*(1-y);
            WzC[1]=(0.75-(0.5-z)*(0.5-z))*effect3D;  
            WzC[2]=0.5*z*z*effect3D;
            WzC[0]=1.0-WzC[1]-WzC[2]; //0.5*(1-z)*(1-z);

            
            i1=(int)(i+x+0.5);
            j1=(int)(j+y+0.5);
            k1=(int)(k+z+0.5);
            x=i+x-i1;
            y=j+y-j1;
            z=k+z-k1;
            Wx[1]=0.75-x*x;
            Wx[2]=0.5*(x+0.5)*(x+0.5);
            Wx[0]=1.0-Wx[1]-Wx[2]; //0.5*(0.5-x)*(0.5-x);
            Wy[1]=(0.75-y*y)*effect2D;
            Wy[2]=0.5*(y+0.5)*(y+0.5)*effect2D;
            Wy[0]=1.0-Wy[1]-Wy[2]; //0.5*(0.5-y)*(0.5-y);
            Wz[1]=(0.75-z*z)*effect3D;
            Wz[2]=0.5*(z+0.5)*(z+0.5)*effect3D;
            Wz[0]=1.0-Wz[1]-Wz[2]; //0.5*(0.5-z)*(0.5-z);

            E2=E3=B2=B3=E1=B1=0.0;
            for(ii=0; ii<iiN; ii++)
              for(jj=0; jj<jjN; jj++)
                for(kk=0; kk<kkN; kk++)
                {
                  E1+=D->Ex[i-1+ii][j1-1+jj+grid2D][k1-1+kk+grid3D]*WxC[ii]*Wy[jj]*Wz[kk];
                  E2+=D->Ey[i1-1+ii][j-1+jj+grid2D][k1-1+kk+grid3D]*Wx[ii]*WyC[jj]*Wz[kk];
                  E3+=D->Ez[i1-1+ii][j1-1+jj+grid2D][k-1+kk+grid3D]*Wx[ii]*Wy[jj]*WzC[kk];
                  B1+=D->Bx[i1-1+ii][j-1+jj+grid2D][k-1+kk+grid3D]*Wx[ii]*WyC[jj]*WzC[kk];
                  B2+=D->By[i-1+ii][j1-1+jj+grid2D][k-1+kk+grid3D]*WxC[ii]*Wy[jj]*WzC[kk];
                  B3+=D->Bz[i-1+ii][j-1+jj+grid2D][k1-1+kk+grid3D]*WxC[ii]*WyC[jj]*Wz[kk];
                }

            p->E1=E1+extE1; p->E2=E2+extE2; p->E3=E3+extE3;
            p->B1=B1+extB1; p->B2=B2+extB2; p->B3=B3+extB3;
      
            p=p->next;
          }
        }		//for(s)  

        //  tracking particle
        s=0;
        p=D->track[i][j][k].head[s]->pt;
        while(p)
        {
          x=p->x;  y=p->y;  z=p->z;
          WxC[2]=0.5*x*x;
          WxC[1]=0.75-(0.5-x)*(0.5-x);
          WxC[0]=1.0-WxC[1]-WxC[2]; //0.5*(1-x)*(1-x);
          WyC[1]=(0.75-(0.5-y)*(0.5-y))*effect2D;
          WyC[2]=0.5*y*y*effect2D;
          WyC[0]=1.0-WyC[1]-WyC[2]; //0.5*(1-y)*(1-y);
          WzC[1]=(0.75-(0.5-z)*(0.5-z))*effect3D;  
          WzC[2]=0.5*z*z*effect3D;
          WzC[0]=1.0-WzC[1]-WzC[2]; //0.5*(1-z)*(1-z);

          
          i1=(int)(i+x+0.5);
          j1=(int)(j+y+0.5);
          k1=(int)(k+z+0.5);
          x=i+x-i1;
          y=j+y-j1;
          z=k+z-k1;
          Wx[1]=0.75-x*x;
          Wx[2]=0.5*(x+0.5)*(x+0.5);
          Wx[0]=1.0-Wx[1]-Wx[2]; //0.5*(0.5-x)*(0.5-x);
          Wy[1]=(0.75-y*y)*effect2D;
          Wy[2]=0.5*(y+0.5)*(y+0.5)*effect2D;
          Wy[0]=1.0-Wy[1]-Wy[2]; //0.5*(0.5-y)*(0.5-y);
          Wz[1]=(0.75-z*z)*effect3D;
          Wz[2]=0.5*(z+0.5)*(z+0.5)*effect3D;
          Wz[0]=1.0-Wz[1]-Wz[2]; //0.5*(0.5-z)*(0.5-z);

          E2=E3=B2=B3=E1=B1=0.0;
          for(ii=0; ii<iiN; ii++)
            for(jj=0; jj<jjN; jj++)
              for(kk=0; kk<kkN; kk++)
              {
                E1+=D->Ex[i-1+ii][j1-1+jj+grid2D][k1-1+kk+grid3D]*WxC[ii]*Wy[jj]*Wz[kk];
                E2+=D->Ey[i1-1+ii][j-1+jj+grid2D][k1-1+kk+grid3D]*Wx[ii]*WyC[jj]*Wz[kk];
                E3+=D->Ez[i1-1+ii][j1-1+jj+grid2D][k-1+kk+grid3D]*Wx[ii]*Wy[jj]*WzC[kk];
                B1+=D->Bx[i1-1+ii][j-1+jj+grid2D][k-1+kk+grid3D]*Wx[ii]*WyC[jj]*WzC[kk];
                B2+=D->By[i-1+ii][j1-1+jj+grid2D][k-1+kk+grid3D]*WxC[ii]*Wy[jj]*WzC[kk];
                B3+=D->Bz[i-1+ii][j-1+jj+grid2D][k1-1+kk+grid3D]*WxC[ii]*WyC[jj]*Wz[kk];
              }

          p->E1=E1+extE1; p->E2=E2+extE2; p->E3=E3+extE3;
          p->B1=B1+extB1; p->B2=B2+extB2; p->B3=B3+extB3;
    
          p=p->next;
        }

      }		   //for(i,j)

}

double secondInterpolation(double ***F,double ***Fxyz,double **Fxy,double *Fx,int i,int j,int k,double x,double y,double z,
                          int jjN,int kkN,double coef2D,double coef3D,int grid2D,int grid3D)
{
  int ii,jj,kk;
  double a,b,c,tmp;

  for(ii=0; ii<3; ii++)
    for(jj=0; jj<jjN; jj++)
      for(kk=0; kk<kkN; kk++) {
        Fxyz[ii][jj][kk]=F[ii+i-1][jj+j-1+grid2D][kk+k-1+grid3D];
      }

  for(ii=0; ii<3; ii++)
    for(jj=0; jj<jjN; jj++) {
      c=Fxyz[ii][jj][0];
      a=0.5*(Fxyz[ii][jj][2]-2*Fxyz[ii][jj][1]+c);
      b=Fxyz[ii][jj][1]-c-a;
      Fxy[ii][jj]=(a*z*z+b*z)*coef3D+c;
      tmp=Fxy[ii][jj];
      if(isnan(tmp)) printf("F=%g,a=%g, b=%g, c=%g,coef3D=%g,z=%g\n",tmp,a,b,c,coef3D,z);
    }

  for(ii=0; ii<3; ii++) {
    c=Fxy[ii][0];
    a=0.5*(Fxy[ii][2]-2*Fxy[ii][1]+c);
    b=Fxy[ii][1]-c-a;
    Fx[ii]=(a*y*y+b*y)*coef2D+c;
  }

  c=Fx[0];
  a=0.5*(Fx[2]-2*Fx[1]+c);
  b=Fx[1]-c-a;
  
//  printf("result=%g,Fx[0]=%g,Fx[1]=%g,Fx[2]=%g\n", a*x*x+b*x+c,Fx[0],Fx[1],Fx[2])  ;

  return a*x*x+b*x+c;
}
