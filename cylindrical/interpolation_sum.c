#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"
#include "math.h"

void interpolation_Split_1st(Domain *D,External *Ext,int iteration);
void interpolation_Split_2nd(Domain *D,External *Ext,int iteration);
void interpolation_Yee_1st(Domain *D,External *Ext);
double interpolation_2nd(double **F,double z,double r);

void interpolation(Domain *D,External *Ext,int iteration)
{
  switch((D->fieldType-1)*2+(D->interpolationType-1)) {  
  case ((Yee-1)*2+(FIRST-1)) :
  case ((NoCherenkov-1)*2+(FIRST-1)) :
    interpolation_Yee_1st(D,Ext);	//2D
    break;
  case ((Split-1)*2+(FIRST-1)) :
    interpolation_Split_1st(D,Ext,iteration);
    break;
  case ((Split-1)*2+(SECOND-1)) :
    interpolation_Split_2nd(D,Ext,iteration);
    break;    
  default :
    printf("In interpolation, what interpolationType(%d)?\n",D->interpolationType);
  }
}

void interpolation_Split_1st(Domain *D,External *Ext,int iteration)
{
  int ii,jj,i,j,m,i1,j1,istart,iend,jstart,jend,s,rank,numMode;
  double coss[D->numMode],sins[D->numMode];
  double coss2[D->numMode],sins2[D->numMode];
  double extE1,extE2,extE3,extB1,extB2,extB3,x,y,z,R,r,invR;
  double wz[2],wr[2],WZ[2],WR[2],y1,y2;
  ptclList *p;
  int myrank, nprocs;      
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  numMode=D->numMode;
	double Pr[2][2],Pl[2][2],Sr[2][2],Sl[2][2],Ez[2][2],Bz[2][2];
  double Pr1[2],Pl1[2],Sr1[2],Sl1[2],Pr2[2],Pl2[2],Sr2[2],Sl2[2];
  double ex,ey,bx,by,ez,er,ep,bz,br,bp,er1,er2,ep1,ep2,br1,br2,bp1,bp2;

  Particle **particle;
  particle=D->particle;

  extE1=Ext->E1;   extE2=Ext->E2;
  extE3=Ext->E3;   extB1=Ext->B1;
  extB2=Ext->B2;   extB3=Ext->B3;

  for(i=istart; i<iend; i++)
    for(j=jstart+1; j<jend; j++) 
      for(s=0; s<D->nSpecies; s++)
      {
        p=particle[i][j].head[s]->pt;
        while(p) {
          z=p->z;  x=p->x; y=p->y;
          R=sqrt(x*x+y*y); invR=1.0/R;
          r=R-(j-jstart);

          i1=((int)(i+z+0.5));
          j1=((int)(j+r+0.5));

          WZ[1]=z+0.5-((int)(z+0.5));  WZ[0]=1.0-WZ[1];
          // WR[1]=(R*R-(j1-jstart-0.5)*(j1-jstart-0.5))/(2.0*(j1-jstart)); WR[0]=1.0-WR[1];
          WR[1]=r+0.5-((int)(r+0.5));  WR[0]=1.0-WR[1];

          wz[1]=z;             wz[0]=1.0-wz[1];
          wr[1]=r;             wr[0]=1.0-wr[1];
          // wr[1]=(R*R-(j-jstart)*(j-jstart))/(2.0*(j-jstart)+1.0);   wr[0]=1.0-wr[1];

          coss[0]=1.0; sins[0]=0.0;
          coss[1]=x*invR; sins[1]=y*invR;
          for(m=2; m<numMode; m++) {
            coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
            sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
          }

          ez=bz=er=br=ep=bp=0.0;
          m=0;
            for(ii=0; ii<2; ii++)
		          for(jj=0; jj<2; jj++) {
			          Pr[ii][jj]=D->PrR[m][i1-1+ii][j1-1+jj];
			          Sr[ii][jj]=D->SrR[m][i1-1+ii][j1-1+jj];
			          Pl[ii][jj]=D->PlR[m][i1-1+ii][j1-1+jj];
			          Sl[ii][jj]=D->SlR[m][i1-1+ii][j1-1+jj];
                Ez[ii][jj]=D->EzNowR[m][i+ii][j+jj];
                Bz[ii][jj]=D->BzNowR[m][i+ii][j+jj];
				      }
          for(m=1; m<numMode; m++) 
		        for(ii=0; ii<2; ii++)
		          for(jj=0; jj<2; jj++) {
			          Pr[ii][jj]+=D->PrR[m][i1-1+ii][j1-1+jj]*coss[m];
			          Pr[ii][jj]-=D->PrI[m][i1-1+ii][j1-1+jj]*sins[m];
			          Sr[ii][jj]+=D->SrR[m][i1-1+ii][j1-1+jj]*coss[m];
			          Sr[ii][jj]-=D->SrI[m][i1-1+ii][j1-1+jj]*sins[m];
			          Pl[ii][jj]+=D->PlR[m][i1-1+ii][j1-1+jj]*coss[m];
			          Pl[ii][jj]-=D->PlI[m][i1-1+ii][j1-1+jj]*sins[m];
			          Sl[ii][jj]+=D->SlR[m][i1-1+ii][j1-1+jj]*coss[m];
			          Sl[ii][jj]-=D->SlI[m][i1-1+ii][j1-1+jj]*sins[m];
                Ez[ii][jj]+=D->EzNowR[m][i+ii][j+jj]*coss[m];
                Ez[ii][jj]-=D->EzNowI[m][i+ii][j+jj]*sins[m];
                Bz[ii][jj]+=D->BzNowR[m][i+ii][j+jj]*coss[m];
                Bz[ii][jj]-=D->BzNowI[m][i+ii][j+jj]*sins[m];
				      }
	     			  
		      for(ii=0; ii<2; ii++)
		        for(jj=0; jj<2; jj++) {
              ez+=wz[ii]*wr[jj]*Ez[ii][jj];
			    	  bz+=wz[ii]*wr[jj]*Bz[ii][jj];
			    	  er+=WZ[ii]*WR[jj]*(Pr[ii][jj]+Pl[ii][jj])*0.5;
			    	  ep+=WZ[ii]*WR[jj]*(Sl[ii][jj]+Sr[ii][jj])*0.5;
			    	  br+=WZ[ii]*WR[jj]*(Sl[ii][jj]-Sr[ii][jj])*0.5;
			    	  bp+=WZ[ii]*WR[jj]*(Pr[ii][jj]-Pl[ii][jj])*0.5;
			      }

          ex=er*coss[1]-ep*sins[1];
          ey=er*sins[1]+ep*coss[1];
          ex=br*coss[1]-bp*sins[1];
          ey=br*sins[1]+bp*coss[1];          

          p->Ez=ez; p->Ex=ex; p->Ey=ey;
          p->Bz=bz; p->Bx=bx; p->By=by;

          p=p->next;
        }
      }		//for(s)        

  j=jstart;
    for(i=istart; i<iend; i++)
      for(s=0; s<D->nSpecies; s++)
      {
        p=particle[i][j].head[s]->pt;
        while(p) {
          z=p->z;  x=p->x; y=p->y;
          r=sqrt(x*x+y*y); invR=1.0/r;

          coss[0]=1.0; sins[0]=0.0;
          coss[1]=x*invR; sins[1]=y*invR;
          for(m=2; m<numMode; m++) {
            coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
            sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
          }
          coss2[0]=1.0; sins2[0]=0.0;
          coss2[1]=-x*invR; sins2[1]=-y*invR;
          for(m=2; m<numMode; m++) {
            coss2[m]=coss2[m-1]*coss2[1]-sins2[m-1]*sins2[1];
            sins2[m]=sins2[m-1]*coss2[1]+coss2[m-1]*sins2[1];
          }

          i1=((int)(i+z+0.5));
          j1=((int)(j+r+0.5));

          wr[1]=r;                     wr[0]=1.0-wr[1];
          // wr[1]=r*r;                   wr[0]=1.0-wr[1];          
          wz[1]=z;                     wz[0]=1.0-wz[1];
          WZ[1]=z+0.5-((int)(z+0.5));  WZ[0]=1.0-WZ[1];
          
          ez=bz=er=br=ep=bp=0.0;
          er1=ep1=er2=ep2=br1=bp1=br2=bp2=0.0;

          //Pr,Pl,Sr,Sl
          if(j1==jstart) {
            m=0;  jj=1;
                for(ii=0; ii<2; ii++) {
			            Pr1[ii]=D->PrR[m][i1-1+ii][j1-1+jj];
			            Sr1[ii]=D->SrR[m][i1-1+ii][j1-1+jj];
			            Pl1[ii]=D->PlR[m][i1-1+ii][j1-1+jj];
			            Sl1[ii]=D->SlR[m][i1-1+ii][j1-1+jj];
                  Pr2[ii]=Pr1[ii];            
                  Sr2[ii]=Sr1[ii];
                  Pl2[ii]=Pl1[ii];            
                  Sl2[ii]=Sl1[ii];
                }
            for(m=1; m<numMode; m++) 
		          for(ii=0; ii<2; ii++) {
			            Pr1[ii]+=D->PrR[m][i1-1+ii][j1-1+jj]*coss[m];
			            Pr1[ii]-=D->PrI[m][i1-1+ii][j1-1+jj]*sins[m];
			            Sr1[ii]+=D->SrR[m][i1-1+ii][j1-1+jj]*coss[m];
			            Sr1[ii]-=D->SrI[m][i1-1+ii][j1-1+jj]*sins[m];
			            Pl1[ii]+=D->PlR[m][i1-1+ii][j1-1+jj]*coss[m];
			            Pl1[ii]-=D->PlI[m][i1-1+ii][j1-1+jj]*sins[m];
			            Sl1[ii]+=D->SlR[m][i1-1+ii][j1-1+jj]*coss[m];
			            Sl1[ii]-=D->SlI[m][i1-1+ii][j1-1+jj]*sins[m];
			            Pr2[ii]+=D->PrR[m][i1-1+ii][j1-1+jj]*coss2[m];
			            Pr2[ii]-=D->PrI[m][i1-1+ii][j1-1+jj]*sins2[m];
			            Sr2[ii]+=D->SrR[m][i1-1+ii][j1-1+jj]*coss2[m];
			            Sr2[ii]-=D->SrI[m][i1-1+ii][j1-1+jj]*sins2[m];
			            Pl2[ii]+=D->PlR[m][i1-1+ii][j1-1+jj]*coss2[m];
			            Pl2[ii]-=D->PlI[m][i1-1+ii][j1-1+jj]*sins2[m];
			            Sl2[ii]+=D->SlR[m][i1-1+ii][j1-1+jj]*coss2[m];
			            Sl2[ii]-=D->SlI[m][i1-1+ii][j1-1+jj]*sins2[m];                  
				        }   
		        for(ii=0; ii<2; ii++) {
			      	  er1+=WZ[ii]*(Pr1[ii]+Pl1[ii])*0.5;
			      	  ep1+=WZ[ii]*(Sl1[ii]+Sr1[ii])*0.5;
			      	  br1+=WZ[ii]*(Sl1[ii]-Sr1[ii])*0.5;
			      	  bp1+=WZ[ii]*(Pr1[ii]-Pl1[ii])*0.5;
			      	  er2+=WZ[ii]*(Pr2[ii]+Pl2[ii])*0.5;
			      	  ep2+=WZ[ii]*(Sl2[ii]+Sr2[ii])*0.5;
			      	  br2+=WZ[ii]*(Sl2[ii]-Sr2[ii])*0.5;
			      	  bp2+=WZ[ii]*(Pr2[ii]-Pl2[ii])*0.5;
			        }                             
            er=(0.5+r)*er1+(0.5-r)*er2;
            ep=(0.5+r)*ep1+(0.5-r)*ep2;
            br=(0.5+r)*br1+(0.5-r)*br2;
            bp=(0.5+r)*bp1+(0.5-r)*bp2;

          } else {
            WR[1]=r+0.5-((int)(r+0.5));  WR[0]=1.0-WR[1];
            // WR[1]=(r*r-(j1-jstart-0.5)*(j1-jstart-0.5))/(2.0*(j1-jstart)); WR[0]=1.0-WR[1];

            m=0;
              for(ii=0; ii<2; ii++)
		            for(jj=0; jj<2; jj++) {
			            Pr[ii][jj]=D->PrR[m][i1-1+ii][j1-1+jj];
			            Sr[ii][jj]=D->SrR[m][i1-1+ii][j1-1+jj];
			            Pl[ii][jj]=D->PlR[m][i1-1+ii][j1-1+jj];
			            Sl[ii][jj]=D->SlR[m][i1-1+ii][j1-1+jj];
				        }
            for(m=1; m<numMode; m++) 
		          for(ii=0; ii<2; ii++)
		            for(jj=0; jj<2; jj++) {
			            Pr[ii][jj]+=D->PrR[m][i1-1+ii][j1-1+jj]*coss[m];
			            Pr[ii][jj]-=D->PrI[m][i1-1+ii][j1-1+jj]*sins[m];
			            Sr[ii][jj]+=D->SrR[m][i1-1+ii][j1-1+jj]*coss[m];
			            Sr[ii][jj]-=D->SrI[m][i1-1+ii][j1-1+jj]*sins[m];
			            Pl[ii][jj]+=D->PlR[m][i1-1+ii][j1-1+jj]*coss[m];
			            Pl[ii][jj]-=D->PlI[m][i1-1+ii][j1-1+jj]*sins[m];
			            Sl[ii][jj]+=D->SlR[m][i1-1+ii][j1-1+jj]*coss[m];
			            Sl[ii][jj]-=D->SlI[m][i1-1+ii][j1-1+jj]*sins[m];
				        }
  
		        for(ii=0; ii<2; ii++)
		          for(jj=0; jj<2; jj++) {
			      	  er+=WZ[ii]*WR[jj]*(Pr[ii][jj]+Pl[ii][jj])*0.5;
			      	  ep+=WZ[ii]*WR[jj]*(Sl[ii][jj]+Sr[ii][jj])*0.5;
			      	  br+=WZ[ii]*WR[jj]*(Sl[ii][jj]-Sr[ii][jj])*0.5;
			      	  bp+=WZ[ii]*WR[jj]*(Pr[ii][jj]-Pl[ii][jj])*0.5;
			        }
			    }
          // Ez, Bz
          m=0;
            for(ii=0; ii<2; ii++)
		          for(jj=0; jj<2; jj++) {
                Ez[ii][jj]=D->EzNowR[m][i+ii][j+jj];
                Bz[ii][jj]=D->BzNowR[m][i+ii][j+jj];
				      }
          for(m=1; m<numMode; m++) 
		        for(ii=0; ii<2; ii++)
		          for(jj=0; jj<2; jj++) {
                Ez[ii][jj]+=D->EzNowR[m][i+ii][j+jj]*coss[m];
                Ez[ii][jj]-=D->EzNowI[m][i+ii][j+jj]*sins[m];
                Bz[ii][jj]+=D->BzNowR[m][i+ii][j+jj]*coss[m];
                Bz[ii][jj]-=D->BzNowI[m][i+ii][j+jj]*sins[m];
				      }	     			  
		      for(ii=0; ii<2; ii++)
		        for(jj=0; jj<2; jj++) {
              ez+=wz[ii]*wr[jj]*Ez[ii][jj];
			    	  bz+=wz[ii]*wr[jj]*Bz[ii][jj];
			      }

          ex=er*coss[1]-ep*sins[1];
          ey=er*sins[1]+ep*coss[1];
          bx=br*coss[1]-bp*sins[1];
          by=br*sins[1]+bp*coss[1];          

          p->Ez=ez; p->Ex=ex; p->Ey=ey;
          p->Bz=bz; p->Bx=bx; p->By=by;
          p=p->next;
        }
      }		//for(s)        

}

void interpolation_Split_2nd(Domain *D,External *Ext,int iteration)
{
  int ii,jj,i,j,m,i1,j1,istart,iend,jstart,jend,s,rank,numMode,intI,intJ;
  double Bz,Br,Bp,Bx,By,Ez,Er,Ep,Ex,Ey,Pr,Pl,Sr,Sl;
  double coss[D->numMode],sins[D->numMode];
  double coss2[D->numMode],sins2[D->numMode];
  double extE1,extE2,extE3,extB1,extB2,extB3,x,y,z,R,r,invR;
  double Pr0,Pr1,Pr2,Pr3,Pl0,Pl1,Pl2,Pl3;
  double Sr0,Sr1,Sr2,Sr3,Sl0,Sl1,Sl2,Sl3;
  double Er1,Er2,Ep1,Ep2,Br1,Br2,Bp1,Bp2,BB,EE,cosTh,sinTh;

   double wz[2],wr[2],WZ[2],WR[2],y1,y2;

  ptclList *p;  
  int myrank, nprocs;    

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  numMode=D->numMode;
	double **PrM,**PlM,**SrM,**SlM,**EzM,**BzM;
	double PrR[numMode][2][2],PlR[numMode][2][2],SrR[numMode][2][2],SlR[numMode][2][2];
	double PrI[numMode][2][2],PlI[numMode][2][2],SrI[numMode][2][2],SlI[numMode][2][2];


  Particle **particle;
  particle=D->particle;

  extE1=Ext->E1;   extE2=Ext->E2;
  extE3=Ext->E3;   extB1=Ext->B1;
  extB2=Ext->B2;   extB3=Ext->B3;

  PrM=(double **)malloc(3*sizeof(double *));
  PlM=(double **)malloc(3*sizeof(double *));  
  SrM=(double **)malloc(3*sizeof(double *));
  SlM=(double **)malloc(3*sizeof(double *));  
  EzM=(double **)malloc(3*sizeof(double *)); 
  BzM=(double **)malloc(3*sizeof(double *)); 
  for(i=0; i<3; i++) {
    PrM[i]=(double *)malloc(3*sizeof(double ));
    PlM[i]=(double *)malloc(3*sizeof(double ));  
    SrM[i]=(double *)malloc(3*sizeof(double ));
    SlM[i]=(double *)malloc(3*sizeof(double ));  
    EzM[i]=(double *)malloc(3*sizeof(double )); 
    BzM[i]=(double *)malloc(3*sizeof(double )); 
  }

  for(i=istart; i<iend; i++)
    for(j=jstart+1; j<jend; j++) 
      for(s=0; s<D->nSpecies; s++)
      {
        p=particle[i][j].head[s]->pt;
        while(p) {
          z=p->z;  x=p->x; y=p->y;
          R=sqrt(x*x+y*y); invR=1.0/R;
          r=R-(j-jstart);

          i1=((int)(i+z+0.5));
          j1=((int)(j+r+0.5));
          intI=(int)(z+0.5);
          intJ=(int)(r+0.5);

          coss[0]=1.0; sins[0]=0.0;
          coss[1]=x*invR; sins[1]=y*invR;
          for(m=2; m<numMode; m++) {
            coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
            sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
          }

          m=0; 
		        for(ii=0; ii<3; ii++)
		          for(jj=0; jj<3; jj++) {
			          PrM[ii][jj]=D->PrR[m][i-1+ii][j-1+jj];          
			          SrM[ii][jj]=D->SrR[m][i-1+ii][j-1+jj];
			          PlM[ii][jj]=D->PlR[m][i-1+ii][j-1+jj];
			          SlM[ii][jj]=D->SlR[m][i-1+ii][j-1+jj];
			          EzM[ii][jj]=D->EzNowR[m][i1-1+ii][j1-1+jj];          
			          BzM[ii][jj]=D->BzNowR[m][i1-1+ii][j1-1+jj];
				      }
          for(m=1; m<numMode; m++) 
		        for(ii=0; ii<3; ii++)
		          for(jj=0; jj<3; jj++) {
			          PrM[ii][jj]+=D->PrR[m][i-1+ii][j-1+jj]*coss[m];
			          PrM[ii][jj]-=D->PrI[m][i-1+ii][j-1+jj]*sins[m];
			          SrM[ii][jj]+=D->SrR[m][i-1+ii][j-1+jj]*coss[m];
			          SrM[ii][jj]-=D->SrI[m][i-1+ii][j-1+jj]*sins[m];
			          PlM[ii][jj]+=D->PlR[m][i-1+ii][j-1+jj]*coss[m];
			          PlM[ii][jj]-=D->PlI[m][i-1+ii][j-1+jj]*sins[m];
			          SlM[ii][jj]+=D->SlR[m][i-1+ii][j-1+jj]*coss[m];
			          SlM[ii][jj]-=D->SlI[m][i-1+ii][j-1+jj]*sins[m];
			          EzM[ii][jj]+=D->EzNowR[m][i1-1+ii][j1-1+jj]*coss[m];
			          EzM[ii][jj]-=D->EzNowI[m][i1-1+ii][j1-1+jj]*sins[m];
			          BzM[ii][jj]+=D->BzNowR[m][i1-1+ii][j1-1+jj]*coss[m];
			          BzM[ii][jj]-=D->BzNowI[m][i1-1+ii][j1-1+jj]*sins[m];
				      }
          Pr=interpolation_2nd(PrM,z-0.5,r-0.5);
          Pl=interpolation_2nd(PlM,z-0.5,r-0.5);
          Sr=interpolation_2nd(SrM,z-0.5,r-0.5);
          Sl=interpolation_2nd(SlM,z-0.5,r-0.5);
          Ez=interpolation_2nd(EzM,z-intI,r-intJ);
          Bz=interpolation_2nd(BzM,z-intI,r-intJ);   

          // Ex=0.5*(Pr+Pl)*coss[1]-0.5*(Sl+Sr)*sins[1];
          // Ey=0.5*(Pr+Pl)*sins[1]+0.5*(Sl+Sr)*coss[1];
          // Bx=0.5*(Sl-Sr)*coss[1]-0.5*(Pr-Pl)*sins[1];
          // By=0.5*(Sl-Sr)*sins[1]+0.5*(Pr-Pl)*coss[1];

          EE=sqrt(Er*Er+Ep*Ep);          
          if(EE==0) { cosTh=0.0; sinTh=0.0; }
          else      { cosTh=Er/EE; sinTh=Ep/EE; } 
          Ex=EE*(coss[1]*cosTh-sins[1]*sinTh);
          Ey=EE*(sins[1]*cosTh+coss[1]*sinTh);
          BB=sqrt(Br*Br+Bp*Bp);          
          if(BB==0) { cosTh=0.0; sinTh=0.0; }
          else      { cosTh=Br/BB; sinTh=Bp/BB; }          
          Bx=BB*(coss[1]*cosTh-sins[1]*sinTh);
          By=BB*(sins[1]*cosTh+coss[1]*sinTh);          

          p->Ez=Ez; p->Ex=Ex; p->Ey=Ey;
          p->Bz=Bz; p->Bx=Bx; p->By=By;

          p=p->next;
        }
      }		//for(s)        

  j=jstart;
    for(i=istart; i<iend; i++)
      for(s=0; s<D->nSpecies; s++)
      {
        p=particle[i][j].head[s]->pt;
        while(p) {
          z=p->z;  x=p->x; y=p->y;
          r=sqrt(x*x+y*y); invR=1.0/r;

          coss[0]=1.0; sins[0]=0.0;
          coss[1]=x*invR; sins[1]=y*invR;
          for(m=2; m<numMode; m++) {
            coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
            sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
          }
          coss2[0]=1.0; sins2[0]=0.0;
          coss2[1]=-x*invR; sins2[1]=-y*invR;
          for(m=2; m<numMode; m++) {
            coss2[m]=coss2[m-1]*coss2[1]-sins2[m-1]*sins2[1];
            sins2[m]=sins2[m-1]*coss2[1]+coss2[m-1]*sins2[1];
          }

          i1=((int)(i+z+0.5));
          j1=((int)(j+r+0.5));

          wr[1]=r;                     wr[0]=1.0-wr[1];
          wz[1]=z;                     wz[0]=1.0-wz[1];
          WZ[1]=z+0.5-((int)(z+0.5));  WZ[0]=1.0-WZ[1];
          
          Ez=Bz=Er=Br=Ep=Bp=0.0;
          for(m=0; m<numMode; m++) 
			      for(ii=0; ii<2; ii++)
			        for(jj=0; jj<2; jj++) {
				        PrR[m][ii][jj]=D->PrR[m][i1-1+ii][j1-1+jj];
				        PrI[m][ii][jj]=D->PrI[m][i1-1+ii][j1-1+jj];
				        SrR[m][ii][jj]=D->SrR[m][i1-1+ii][j1-1+jj];
				        SrI[m][ii][jj]=D->SrI[m][i1-1+ii][j1-1+jj];
				        PlR[m][ii][jj]=D->PlR[m][i1-1+ii][j1-1+jj];
				        PlI[m][ii][jj]=D->PlI[m][i1-1+ii][j1-1+jj];
				        SlR[m][ii][jj]=D->SlR[m][i1-1+ii][j1-1+jj];
				        SlI[m][ii][jj]=D->SlI[m][i1-1+ii][j1-1+jj];
				      }

          //Pr,Pl,Sr,Sl
          if(j1==jstart) {
            WR[1]=r/0.5;     WR[0]=1.0-WR[1];

            // m=0; jj=1;
			      // for(ii=0; ii<2; ii++) {
				    //   Er+=WZ[ii]*WR[jj]*(PrR[m][ii][jj]+PlR[m][ii][jj])*0.5;
				    //   Ep+=WZ[ii]*WR[jj]*(SlR[m][ii][jj]+SrR[m][ii][jj])*0.5;
				    //   Br+=WZ[ii]*WR[jj]*(SlR[m][ii][jj]-SrR[m][ii][jj])*0.5;
				    //   Bp+=WZ[ii]*WR[jj]*(PrR[m][ii][jj]-PlR[m][ii][jj])*0.5;
				    // }
            // m=1;
			      // for(ii=0; ii<2; ii++) {
  				  //   Er+=WZ[ii]*WR[jj]*(PrR[m][ii][jj]+PlR[m][ii][jj])*0.5*coss[m];
  				  //   Er-=WZ[ii]*WR[jj]*(PrI[m][ii][jj]+PlI[m][ii][jj])*0.5*sins[m];
				    //   Ep+=WZ[ii]*WR[jj]*(SlR[m][ii][jj]+SrR[m][ii][jj])*0.5*coss[m];
				    //   Ep-=WZ[ii]*WR[jj]*(SlI[m][ii][jj]+SrI[m][ii][jj])*0.5*sins[m];
				    //   Br+=WZ[ii]*WR[jj]*(SlR[m][ii][jj]-SrR[m][ii][jj])*0.5*coss[m];
				    //   Br-=WZ[ii]*WR[jj]*(SlI[m][ii][jj]-SrI[m][ii][jj])*0.5*sins[m];
				    //   Bp+=WZ[ii]*WR[jj]*(PrR[m][ii][jj]-PlR[m][ii][jj])*0.5*coss[m];
				    //   Bp-=WZ[ii]*WR[jj]*(PrI[m][ii][jj]-PlI[m][ii][jj])*0.5*sins[m];
				    // }
            // for(m=2; m<numMode; m++) 
  			    //   for(ii=0; ii<2; ii++) {
  				  //     Er+=WZ[ii]*WR[jj]*(PrR[m][ii][jj]+PlR[m][ii][jj])*0.5*coss[m];
  				  //     Er-=WZ[ii]*WR[jj]*(PrI[m][ii][jj]+PlI[m][ii][jj])*0.5*sins[m];
				    //     Ep+=WZ[ii]*WR[jj]*(SlR[m][ii][jj]+SrR[m][ii][jj])*0.5*coss[m];
        		//     Ep-=WZ[ii]*WR[jj]*(SlI[m][ii][jj]+SrI[m][ii][jj])*0.5*sins[m];
				    //     Br+=WZ[ii]*WR[jj]*(SlR[m][ii][jj]-SrR[m][ii][jj])*0.5*coss[m];
				    //     Br-=WZ[ii]*WR[jj]*(SlI[m][ii][jj]-SrI[m][ii][jj])*0.5*sins[m];
				    //     Bp+=WZ[ii]*WR[jj]*(PrR[m][ii][jj]-PlR[m][ii][jj])*0.5*coss[m];
				    //     Bp-=WZ[ii]*WR[jj]*(PrI[m][ii][jj]-PlI[m][ii][jj])*0.5*sins[m];
            //   }

			      Pr1=Pl1=Sr1=Sl1=0.0;
			      Pr2=Pl2=Sr2=Sl2=0.0;
            m=0; jj=1;
			      for(ii=0; ii<2; ii++) {
				  	  Pr1+=WZ[ii]*D->PrR[m][i1-1+ii][j1-1+jj];
				  	  Pl1+=WZ[ii]*D->PlR[m][i1-1+ii][j1-1+jj];
				  	  Sr1+=WZ[ii]*D->SrR[m][i1-1+ii][j1-1+jj];
				  	  Sl1+=WZ[ii]*D->SlR[m][i1-1+ii][j1-1+jj];
				    }
            Pr2=Pr1; Pl2=Pl2; Sr2=Sr1; Sl2=Sl1;
            for(m=1; m<numMode; m++) 
  			      for(ii=0; ii<2; ii++) { 
  			  	    Pr1+=WZ[ii]*D->PrR[m][i1-1+ii][j1-1+jj]*coss[m];
  			  	    Pr1-=WZ[ii]*D->PrI[m][i1-1+ii][j1-1+jj]*sins[m];
				  	    Pl1+=WZ[ii]*D->PlR[m][i1-1+ii][j1-1+jj]*coss[m];
				  	    Pl1-=WZ[ii]*D->PlI[m][i1-1+ii][j1-1+jj]*sins[m];
				  	    Sr1+=WZ[ii]*D->SrR[m][i1-1+ii][j1-1+jj]*coss[m];
				  	    Sr1-=WZ[ii]*D->SrI[m][i1-1+ii][j1-1+jj]*sins[m];
				  	    Sl1+=WZ[ii]*D->SlR[m][i1-1+ii][j1-1+jj]*coss[m];
				  	    Sl1-=WZ[ii]*D->SlI[m][i1-1+ii][j1-1+jj]*sins[m];

  			  	    Pr2+=WZ[ii]*D->PrR[m][i1-1+ii][j1-1+jj]*coss2[m];
  			  	    Pr2-=WZ[ii]*D->PrI[m][i1-1+ii][j1-1+jj]*sins2[m];
				  	    Pl2+=WZ[ii]*D->PlR[m][i1-1+ii][j1-1+jj]*coss2[m];
				  	    Pl2-=WZ[ii]*D->PlI[m][i1-1+ii][j1-1+jj]*sins2[m];
				  	    Sr2+=WZ[ii]*D->SrR[m][i1-1+ii][j1-1+jj]*coss2[m];
				  	    Sr2-=WZ[ii]*D->SrI[m][i1-1+ii][j1-1+jj]*sins2[m];
				  	    Sl2+=WZ[ii]*D->SlR[m][i1-1+ii][j1-1+jj]*coss2[m];
				  	    Sl2-=WZ[ii]*D->SlI[m][i1-1+ii][j1-1+jj]*sins2[m];
              }
            Er=(0.5+r)*0.5*(Pr1+Pl1)+(0.5-r)*0.5*(Pr2+Pl2);
            Ep=(0.5+r)*0.5*(Sr1+Sl1)+(0.5-r)*0.5*(Sr2+Sl2);
            Br=(0.5+r)*0.5*(Sl1-Sr1)+(0.5-r)*0.5*(Sl2-Sr2);
            Bp=(0.5+r)*0.5*(Pr1-Pl1)+(0.5-r)*0.5*(Pr2-Pl2);

          } else {
            WR[1]=r+0.5-((int)(r+0.5));  WR[0]=1.0-WR[1];

            m=0;
			        for(ii=0; ii<2; ii++)
			          for(jj=0; jj<2; jj++) {
					        Er+=WZ[ii]*WR[jj]*(PrR[m][ii][jj]+PlR[m][ii][jj])*0.5;
					        Ep+=WZ[ii]*WR[jj]*(SlR[m][ii][jj]+SrR[m][ii][jj])*0.5;
					        Br+=WZ[ii]*WR[jj]*(SlR[m][ii][jj]-SrR[m][ii][jj])*0.5;
					        Bp+=WZ[ii]*WR[jj]*(PrR[m][ii][jj]-PlR[m][ii][jj])*0.5;
				        }
            for(m=1; m<numMode; m++) 
			        for(ii=0; ii<2; ii++)
			          for(jj=0; jj<2; jj++) {
  					      Er+=WZ[ii]*WR[jj]*(PrR[m][ii][jj]+PlR[m][ii][jj])*0.5*coss[m];
  					      Er-=WZ[ii]*WR[jj]*(PrI[m][ii][jj]+PlI[m][ii][jj])*0.5*sins[m];
					        Ep+=WZ[ii]*WR[jj]*(SlR[m][ii][jj]+SrR[m][ii][jj])*0.5*coss[m];
					        Ep-=WZ[ii]*WR[jj]*(SlI[m][ii][jj]+SrI[m][ii][jj])*0.5*sins[m];
					        Br+=WZ[ii]*WR[jj]*(SlR[m][ii][jj]-SrR[m][ii][jj])*0.5*coss[m];
					        Br-=WZ[ii]*WR[jj]*(SlI[m][ii][jj]-SrI[m][ii][jj])*0.5*sins[m];
					        Bp+=WZ[ii]*WR[jj]*(PrR[m][ii][jj]-PlR[m][ii][jj])*0.5*coss[m];
					        Bp-=WZ[ii]*WR[jj]*(PrI[m][ii][jj]-PlI[m][ii][jj])*0.5*sins[m];
                }
			    }

          m=0;
          for(ii=0; ii<2; ii++)
			      for(jj=0; jj<2; jj++) {
              Ez+=wz[ii]*wr[jj]*D->EzNowR[m][i+ii][j+jj];
					    Bz+=wz[ii]*wr[jj]*D->BzNowR[m][i+ii][j+jj];
				    }
			    jj=1;
          for(m=1; m<numMode; m++) 
            for(ii=0; ii<2; ii++)  {
              Ez+=wz[ii]*wr[jj]*D->EzNowR[m][i+ii][j+jj]*coss[m];
              Ez-=wz[ii]*wr[jj]*D->EzNowI[m][i+ii][j+jj]*sins[m];
              Bz+=wz[ii]*wr[jj]*D->BzNowR[m][i+ii][j+jj]*coss[m];
              Bz-=wz[ii]*wr[jj]*D->BzNowI[m][i+ii][j+jj]*sins[m];
			      }

          // Ex=Er*coss[1]-Ep*sins[1];
          // Ey=Er*sins[1]+Ep*coss[1];
          // Bx=Br*coss[1]-Bp*sins[1];
          // By=Br*sins[1]+Bp*coss[1];

          EE=sqrt(Er*Er+Ep*Ep);          
          if(EE==0) { cosTh=0.0; sinTh=0.0; }
          else      { cosTh=Er/EE; sinTh=Ep/EE; } 
          Ex=EE*(coss[1]*cosTh-sins[1]*sinTh);
          Ey=EE*(sins[1]*cosTh+coss[1]*sinTh);
          BB=sqrt(Br*Br+Bp*Bp);          
          if(BB==0) { cosTh=0.0; sinTh=0.0; }
          else      { cosTh=Br/BB; sinTh=Bp/BB; }          
          Bx=BB*(coss[1]*cosTh-sins[1]*sinTh);
          By=BB*(sins[1]*cosTh+coss[1]*sinTh);   

          p->Ez=Ez; p->Ex=Ex; p->Ey=Ey;
          p->Bz=Bz; p->Bx=Bx; p->By=By;
          p=p->next;
        }
      }		//for(s)    

/*
  j=jstart;
    for(i=istart; i<iend; i++)
      for(s=0; s<D->nSpecies; s++)
      {
        p=particle[i][j].head[s]->pt;
        while(p) {
          z=p->z;  x=p->x; y=p->y;
          r=sqrt(x*x+y*y); invR=1.0/r;

          coss[0]=1.0; sins[0]=0.0;
          coss[1]=x*invR; sins[1]=y*invR;
          for(m=2; m<numMode; m++) {
            coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
            sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
          }
          coss2[0]=1.0; sins2[0]=0.0;
          coss2[1]=-x*invR; sins2[1]=-y*invR;
          for(m=2; m<numMode; m++) {
            coss2[m]=coss2[m-1]*coss2[1]-sins2[m-1]*sins2[1];
            sins2[m]=sins2[m-1]*coss2[1]+coss2[m-1]*sins2[1];
          }

          i1=((int)(i+z+0.5));
          j1=((int)(j+r+0.5));
          intI=(int)(z+0.5);
          intJ=(int)(r+0.5);

          m=0; 
		        for(ii=0; ii<3; ii++) {
              jj=0;
			          PrM[ii][jj]=D->PrR[m][i-1+ii][j+jj];          
			          SrM[ii][jj]=D->SrR[m][i-1+ii][j+jj];
			          PlM[ii][jj]=D->PlR[m][i-1+ii][j+jj];
			          SlM[ii][jj]=D->SlR[m][i-1+ii][j+jj];
		          for(jj=1; jj<3; jj++) {
			          PrM[ii][jj]=D->PrR[m][i-1+ii][j-1+jj];          
			          SrM[ii][jj]=D->SrR[m][i-1+ii][j-1+jj];
			          PlM[ii][jj]=D->PlR[m][i-1+ii][j-1+jj];
			          SlM[ii][jj]=D->SlR[m][i-1+ii][j-1+jj];
				      }
            }
          for(m=1; m<numMode; m++) 
		        for(ii=0; ii<3; ii++) {
              jj=0;
			          PrM[ii][jj]+=D->PrR[m][i-1+ii][j+jj]*coss2[m];
			          PrM[ii][jj]-=D->PrI[m][i-1+ii][j+jj]*sins2[m];
			          SrM[ii][jj]+=D->SrR[m][i-1+ii][j+jj]*coss2[m];
			          SrM[ii][jj]-=D->SrI[m][i-1+ii][j+jj]*sins2[m];
			          PlM[ii][jj]+=D->PlR[m][i-1+ii][j+jj]*coss2[m];
			          PlM[ii][jj]-=D->PlI[m][i-1+ii][j+jj]*sins2[m];
			          SlM[ii][jj]+=D->SlR[m][i-1+ii][j+jj]*coss2[m];
			          SlM[ii][jj]-=D->SlI[m][i-1+ii][j+jj]*sins2[m];
		          for(jj=1; jj<3; jj++) {
			          PrM[ii][jj]+=D->PrR[m][i-1+ii][j-1+jj]*coss[m];
			          PrM[ii][jj]-=D->PrI[m][i-1+ii][j-1+jj]*sins[m];
			          SrM[ii][jj]+=D->SrR[m][i-1+ii][j-1+jj]*coss[m];
			          SrM[ii][jj]-=D->SrI[m][i-1+ii][j-1+jj]*sins[m];
			          PlM[ii][jj]+=D->PlR[m][i-1+ii][j-1+jj]*coss[m];
			          PlM[ii][jj]-=D->PlI[m][i-1+ii][j-1+jj]*sins[m];
			          SlM[ii][jj]+=D->SlR[m][i-1+ii][j-1+jj]*coss[m];
			          SlM[ii][jj]-=D->SlI[m][i-1+ii][j-1+jj]*sins[m];
				      }
            }

          if(j1==jstart) {
            m=0; 
		          for(ii=0; ii<3; ii++) {
		            jj=0;
			            EzM[ii][jj]=D->EzNowR[m][i1-1+ii][j1+1+jj];          
			            BzM[ii][jj]=D->BzNowR[m][i1-1+ii][j1+1+jj];
		            for(jj=1; jj<3; jj++) {
			            EzM[ii][jj]=D->EzNowR[m][i1-1+ii][j1-1+jj];          
			            BzM[ii][jj]=D->BzNowR[m][i1-1+ii][j1-1+jj];
				        }
              }
            for(m=1; m<numMode; m++) 
		          for(ii=0; ii<3; ii++) {
                jj=0;
			            EzM[ii][jj]+=D->EzNowR[m][i1-1+ii][j1+1+jj]*coss2[m];
			            EzM[ii][jj]-=D->EzNowI[m][i1-1+ii][j1+1+jj]*sins2[m];
			            BzM[ii][jj]+=D->BzNowR[m][i1-1+ii][j1+1+jj]*coss2[m];
			            BzM[ii][jj]-=D->BzNowI[m][i1-1+ii][j1+1+jj]*sins2[m];
		            for(jj=1; jj<3; jj++) {
			            EzM[ii][jj]+=D->EzNowR[m][i1-1+ii][j1-1+jj]*coss[m];
			            EzM[ii][jj]-=D->EzNowI[m][i1-1+ii][j1-1+jj]*sins[m];
			            BzM[ii][jj]+=D->BzNowR[m][i1-1+ii][j1-1+jj]*coss[m];
			            BzM[ii][jj]-=D->BzNowI[m][i1-1+ii][j1-1+jj]*sins[m];
				        }
              }
          } else {
            m=0; 
		          for(ii=0; ii<3; ii++)
		            for(jj=0; jj<3; jj++) {
			            EzM[ii][jj]=D->EzNowR[m][i1-1+ii][j1-1+jj];          
			            BzM[ii][jj]=D->BzNowR[m][i1-1+ii][j1-1+jj];
				        }
            for(m=1; m<numMode; m++) 
		          for(ii=0; ii<3; ii++)
		            for(jj=0; jj<3; jj++) {
			            EzM[ii][jj]+=D->EzNowR[m][i1-1+ii][j1-1+jj]*coss[m];
			            EzM[ii][jj]-=D->EzNowI[m][i1-1+ii][j1-1+jj]*sins[m];
			            BzM[ii][jj]+=D->BzNowR[m][i1-1+ii][j1-1+jj]*coss[m];
			            BzM[ii][jj]-=D->BzNowI[m][i1-1+ii][j1-1+jj]*sins[m];
				        }
          }
          Pr=interpolation_2nd(PrM,z-0.5,r-0.5);
          Pl=interpolation_2nd(PlM,z-0.5,r-0.5);
          Sr=interpolation_2nd(SrM,z-0.5,r-0.5);
          Sl=interpolation_2nd(SlM,z-0.5,r-0.5);          
          Ez=interpolation_2nd(EzM,z-intI,r-intJ);
          Bz=interpolation_2nd(BzM,z-intI,r-intJ);          

          Ex=0.5*(Pr+Pl)*coss[1]-0.5*(Sl+Sr)*sins[1];
          Ey=0.5*(Pr+Pl)*sins[1]+0.5*(Sl+Sr)*coss[1];
          Bx=0.5*(Sl-Sr)*coss[1]-0.5*(Pr-Pl)*sins[1];
          By=0.5*(Sl-Sr)*sins[1]+0.5*(Pr-Pl)*coss[1];
          p->Ez=Ez; p->Ex=Ex; p->Ey=Ey;
          p->Bz=Bz; p->Bx=Bx; p->By=By;

          p=p->next;
        }
      }		//for(s)        
*/
  for(i=0; i<3; i++) {
    free(PrM[i]);
    free(PlM[i]);
    free(SrM[i]);
    free(SlM[i]);
    free(EzM[i]);
    free(BzM[i]);
  }
  free(PrM);
  free(PlM);
  free(SrM);
  free(SlM);
  free(EzM);
  free(BzM);  
}


void interpolation_Yee_1st(Domain *D,External *Ext)
{
   int i,j,m,i1,j1,ii,jj,istart,iend,jstart,jend,s,rank,numMode;
   double Bz,Br,Bp,Bx,By,Ez,Er,Ep,Ex,Ey,yR1,yR2,yI1,yI2;
   double coss[D->numMode],sins[D->numMode];
   double extE1,extE2,extE3,extB1,extB2,extB3,x,y,z,R,r,invR;
   double wz[2],wr[2],WZ[2],WR[2];
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   numMode=D->numMode;

   Particle **particle;
   particle=D->particle;

   extE1=Ext->E1;   extE2=Ext->E2;
   extE3=Ext->E3;   extB1=Ext->B1;
   extB2=Ext->B2;   extB3=Ext->B3;

   for(i=istart; i<iend; i++)
     for(j=jstart+1; j<jend; j++) 
     { 
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j].head[s]->pt;
         while(p) {
           z=p->z;  x=p->x; y=p->y;
           R=sqrt(x*x+y*y); invR=1.0/R;
           r=R-(j-jstart);

           wz[1]=z;           wz[0]=1.0-wz[1];
           wr[1]=r;           wr[0]=1.0-wr[1];
           i1=((int)(i+z+0.5));
           j1=((int)(j+r+0.5));
           WZ[1]=z+0.5-((int)(z+0.5));  WZ[0]=1.0-WZ[1];
           WR[1]=r+0.5-((int)(r+0.5));  WR[0]=1.0-WR[1];

           coss[0]=1.0; sins[0]=0.0;
           coss[1]=x*invR; sins[1]=y*invR;
           for(m=2; m<numMode; m++) {
             coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
             sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
           }

           Bz=Br=Bp=Ez=Er=Ep=0.0;
           for(ii=0; ii<2; ii++)
             for(jj=0; jj<2; jj++) {
               Bz+=wz[ii]*WR[jj]*D->BzNowR[0][i+ii][j1-1+jj];
               Br+=WZ[ii]*wr[jj]*D->BrNowR[0][i1-1+ii][j+jj];
               Bp+=WZ[ii]*WR[jj]*D->BpNowR[0][i1-1+ii][j1-1+jj];
               Ez+=WZ[ii]*wr[jj]*D->EzR[0][i1-1+ii][j+jj];
               Er+=wz[ii]*WR[jj]*D->ErR[0][i+ii][j1-1+jj];
               Ep+=wz[ii]*wr[jj]*D->EpR[0][i+ii][j+jj];
             }

           for(m=1; m<numMode; m++) {
             for(ii=0; ii<2; ii++)
               for(jj=0; jj<2; jj++) {
                 Bz+=wz[ii]*WR[jj]*D->BzNowR[m][i+ii][j1-1+jj]*coss[m];
                 Bz-=wz[ii]*WR[jj]*D->BzNowI[m][i+ii][j1-1+jj]*sins[m];
                 Br+=WZ[ii]*wr[jj]*D->BrNowR[m][i1-1+ii][j+jj]*coss[m];
                 Br-=WZ[ii]*wr[jj]*D->BrNowI[m][i1-1+ii][j+jj]*sins[m];
                 Bp+=WZ[ii]*WR[jj]*D->BpNowR[m][i1-1+ii][j1-1+jj]*coss[m];
                 Bp-=WZ[ii]*WR[jj]*D->BpNowI[m][i1-1+ii][j1-1+jj]*sins[m];
                 Ez+=WZ[ii]*wr[jj]*D->EzR[m][i1-1+ii][j+jj]*coss[m];
                 Ez-=WZ[ii]*wr[jj]*D->EzI[m][i1-1+ii][j+jj]*sins[m];
                 Er+=wz[ii]*WR[jj]*D->ErR[m][i+ii][j1-1+jj]*coss[m];
                 Er-=wz[ii]*WR[jj]*D->ErI[m][i+ii][j1-1+jj]*sins[m];
                 Ep+=wz[ii]*wr[jj]*D->EpR[m][i+ii][j+jj]*coss[m];
                 Ep-=wz[ii]*wr[jj]*D->EpI[m][i+ii][j+jj]*sins[m];
               }
           }
           Ex=Er*coss[1]-Ep*sins[1];
           Ey=Er*sins[1]+Ep*coss[1];
           Bx=Br*coss[1]-Bp*sins[1];
           By=Br*sins[1]+Bp*coss[1];

           p->Ez=Ez; p->Ex=Ex; p->Ey=Ey;
           p->Bz=Bz; p->Bx=Bx; p->By=By;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)

   j=jstart;
   for(i=istart; i<iend; i++)
     { 
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j].head[s]->pt;
         while(p) {
           z=p->z;  x=p->x; y=p->y;
           R=sqrt(x*x+y*y); invR=1.0/R;
           r=R-(j-jstart);

           coss[0]=1.0; sins[0]=0.0;
           coss[1]=x*invR; sins[1]=y*invR;
           for(m=2; m<numMode; m++) {
             coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
             sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
           }

           i1=((int)(i+z+0.5));
           j1=((int)(j+r+0.5));
           wz[1]=z;                     wz[0]=1.0-wz[1];
           WZ[1]=z+0.5-((int)(z+0.5));  WZ[0]=1.0-WZ[1];

           Bz=Br=Bp=Ez=Er=Ep=0.0;
           //Bz, Bp, Er
           if(j1==jstart) {
             WR[1]=R*2.0;                 WR[0]=1.0-WR[1];
             m=0;
             for(ii=0; ii<2; ii++) {
               Bz+=wz[ii]*D->BzNowR[m][i+ii][j1];
               Bp+=WZ[ii]*WR[1]*D->BpNowR[m][i1-1+ii][j1];
               Er+=wz[ii]*WR[1]*D->ErR[m][i+ii][j1];
             }
             m=1;
             for(ii=0; ii<2; ii++) {
               Bz+=wz[ii]*WR[1]*D->BzNowR[m][i+ii][j1]*coss[m]
                  -wz[ii]*WR[1]*D->BzNowI[m][i+ii][j1]*sins[m];
               Bp+=WZ[ii]*D->BpNowR[m][i1-1+ii][j1]*coss[m]
                  -WZ[ii]*D->BpNowI[m][i1-1+ii][j1]*sins[m];
               Er+=wz[ii]*D->ErR[m][i+ii][j1]*coss[m]
                  -wz[ii]*D->ErI[m][i+ii][j1]*sins[m];
             }
             for(m=2; m<numMode; m++) 
               for(ii=0; ii<2; ii++) {
                 Bz+=wz[ii]*WR[1]*D->BzNowR[m][i+ii][j1]*coss[m]
                    -wz[ii]*WR[1]*D->BzNowI[m][i+ii][j1]*sins[m];
                 Bp+=WZ[ii]*WR[1]*D->BpNowR[m][i1-1+ii][j1]*coss[m]
                    -WZ[ii]*WR[1]*D->BpNowI[m][i1-1+ii][j1]*sins[m];
                 Er+=wz[ii]*WR[1]*D->ErR[m][i+ii][j1]*coss[m]
                    -wz[ii]*WR[1]*D->ErI[m][i+ii][j1]*sins[m];
               }
           } else {
             WR[1]=r+0.5-((int)(r+0.5));  WR[0]=1.0-WR[1];

             for(ii=0; ii<2; ii++)
               for(jj=0; jj<2; jj++) {
                 Bz+=wz[ii]*WR[jj]*D->BzNowR[0][i+ii][j1-1+jj];
                 Bp+=WZ[ii]*WR[jj]*D->BpNowR[0][i1-1+ii][j1-1+jj];
                 Er+=wz[ii]*WR[jj]*D->ErR[0][i+ii][j1-1+jj];
               }
             for(m=1; m<numMode; m++) 
               for(ii=0; ii<2; ii++)
                 for(jj=0; jj<2; jj++) {
                   Bz+=wz[ii]*WR[jj]*D->BzNowR[m][i+ii][j1-1+jj]*coss[m];
                   Bz-=wz[ii]*WR[jj]*D->BzNowI[m][i+ii][j1-1+jj]*sins[m];
                   Bp+=WZ[ii]*WR[jj]*D->BpNowR[m][i1-1+ii][j1-1+jj]*coss[m];
                   Bp-=WZ[ii]*WR[jj]*D->BpNowI[m][i1-1+ii][j1-1+jj]*sins[m];
                   Er+=wz[ii]*WR[jj]*D->ErR[m][i+ii][j1-1+jj]*coss[m];
                   Er-=wz[ii]*WR[jj]*D->ErI[m][i+ii][j1-1+jj]*sins[m];
                 }             
           }

           //Ez, Ep, Br
           wr[1]=r;                  wr[0]=1.0-wr[1];

           for(ii=0; ii<2; ii++)
             for(jj=0; jj<2; jj++) {
               Br+=WZ[ii]*wr[jj]*D->BrNowR[0][i1-1+ii][j+jj];
               Ez+=WZ[ii]*wr[jj]*D->EzR[0][i1-1+ii][j+jj];
               Ep+=wz[ii]*wr[jj]*D->EpR[0][i+ii][j+jj];
             }

           for(m=1; m<numMode; m++) 
             for(ii=0; ii<2; ii++)
               for(jj=0; jj<2; jj++) {
                 Br+=WZ[ii]*wr[jj]*D->BrNowR[m][i1-1+ii][j+jj]*coss[m]
                    -WZ[ii]*wr[jj]*D->BrNowI[m][i1-1+ii][j+jj]*sins[m];
                 Ez+=WZ[ii]*wr[jj]*D->EzR[m][i1-1+ii][j+jj]*coss[m]
                    -WZ[ii]*wr[jj]*D->EzI[m][i1-1+ii][j+jj]*sins[m];
                 Ep+=wz[ii]*wr[jj]*D->EpR[m][i+ii][j+jj]*coss[m]
                    -wz[ii]*wr[jj]*D->EpI[m][i+ii][j+jj]*sins[m];
               }
           
           Ex=Er*coss[1]-Ep*sins[1];
           Ey=Er*sins[1]+Ep*coss[1];
           Bx=Br*coss[1]-Bp*sins[1];
           By=Br*sins[1]+Bp*coss[1];

           p->Ez=Ez; p->Ex=Ex; p->Ey=Ey;
           p->Bz=Bz; p->Bx=Bx; p->By=By;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)

}

double interpolation_2nd(double **F,double z,double r)
{
  int i,j;
  double Fx[3],a,b,c;

  for(i=0; i<3; i++) {
    c=F[i][1];
    b=0.5*(F[i][2]-F[i][0]);
    a=0.5*(F[i][2]+F[i][0])-c;
    Fx[i]=a*r*r+b*r+c;
  }
  c=Fx[1];
  b=0.5*(Fx[2]-Fx[0]);
  a=0.5*(Fx[2]+Fx[0])-c;

  return a*z*z+b*z+c;
}



/*
void interpolation2D_Yee_Pukhov_2nd(Domain *D,External *Ext)
{
   int s,i,j,k,ii,jj,i,j+1,k1,istart,iend,jstart,jend,kstart,kend;
   double E1,E2,E3,B1,B2,B3,extE1,extE2,extE3,extB1,extB2,extB3;
   double x,y,z,xx,yy,zz,x1,y1,z1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;   extE2=Ext->E2;   extE3=Ext->E3;
   extB1=Ext->B1;   extB2=Ext->B2;   extB3=Ext->B3;
   
   k=k1=0;
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x;  y=p->y; // z=p->z;
           //edge
           ii=(int)(i+x+0.5);           jj=(int)(j+y+0.5);
           xx=1.0+x-((int)(0.5+x));     yy=1.0+y-((int)(0.5+y));
           //side
           i=i;           j+1=j;
           x1=0.5+x;       y1=0.5+y;

           E1=calBi+1D(D->Ex,i,jj,x1,yy);
           E2=calBi+1D(D->Ey,ii,j+1,xx,y1);
           E3=calBi+1D(D->Ez,ii,jj,xx,yy);
           B1=calBi+1D(D->BxNow,ii,j+1,xx,y1);
           B2=calBi+1D(D->ByNow,i,jj,x1,yy);
           B3=calBi+1D(D->BzNow,i,j+1,x1,y1);

           p->E1=E1+extE1; p->E2=E2+extE2; p->E3=E3+extE3;
           p->B1=B1+extB1; p->B2=B2+extB2; p->B3=B3+extB3;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)
}

void interpolation1D_Split_1st(Domain *D,External *Ext)
{
   int i,j,k,i,istart,iend,s,cnt;
   double E1,Pr,Pl,B1,Sr,Sl,extE1,extE2,extE3,extB1,extB2,extB3,x,x1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart; iend=D->iend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;   extE2=Ext->E2;   extE3=Ext->E3;
   extB1=Ext->B1;   extB2=Ext->B2;   extB3=Ext->B3;
   
   j=k=0;
   for(i=istart; i<iend; i++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x; 
           i=(int)(i+x+0.5);
           x1=x+0.5-((int)(x+0.5));

           B1=0;
           E1=(1-x1)*D->Ex[i-1][j][k] + x1*D->Ex[i][j][k];
           Pr=(1-x1)*D->Pr[i-1][j][k] + x1*D->Pr[i][j][k];
           Pl=(1-x1)*D->Pl[i-1][j][k] + x1*D->Pl[i][j][k];
           Sr=(1-x1)*D->Sr[i-1][j][k] + x1*D->Sr[i][j][k];
           Sl=(1-x1)*D->Sl[i-1][j][k] + x1*D->Sl[i][j][k];

           p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
           p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)
}

void interpolation2D_Split_1st(Domain *D,External *Ext)
{
   int i,j,k,i,j+1,k1,istart,iend,jstart,jend,s,cnt;
   double E1,Pr,Pl,B1,Sr,Sl,x,y,z,x1,y1,z1;
   double extE1,extE2,extE3,extB1,extB2,extB3;
   ptclList *p;
   int myrank, nprocs;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;   extE2=Ext->E2;   extE3=Ext->E3;
   extB1=Ext->B1;   extB2=Ext->B2;   extB3=Ext->B3;

   k=k1=0;
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         cnt=0;
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x;  y=p->y; // z=p->z;
           i=(int)(i+x+0.5);           j+1=(int)(j+y+0.5);
           x1=x+0.5-((int)(x+0.5));     y1=y+0.5-((int)(y+0.5));

           B1=(1-x1)*(1-y1)*D->Bx[i-1][j][k1]
             +    x1*(1-y1)*D->Bx[i][j][k1]
             +(1-x1)*    y1*D->Bx[i-1][j+1][k1]
             +    x1*    y1*D->Bx[i][j+1][k1];
           E1=(1-x1)*(1-y)*D->Ex[i-1][j][k]
             +    x1*(1-y)*D->Ex[i][j][k]
             +(1-x1)*    y*D->Ex[i-1][j+1][k]
             +    x1*    y*D->Ex[i][j+1][k];
           Pr=(1-x1)*(1-y1)*D->Pr[i-1][j][k]
             +    x1*(1-y1)*D->Pr[i][j][k]
             +(1-x1)*    y1*D->Pr[i-1][j+1][k]
             +    x1*    y1*D->Pr[i][j+1][k];
           Pl=(1-x1)*(1-y1)*D->Pl[i-1][j][k]
             +    x1*(1-y1)*D->Pl[i][j][k]
             +(1-x1)*    y1*D->Pl[i-1][j+1][k]
             +    x1*    y1*D->Pl[i][j+1][k];
           Sr=(1-x1)*(1-y)*D->Sr[i-1][j][k1]
             +    x1*(1-y)*D->Sr[i][j][k1]
             +(1-x1)*    y*D->Sr[i-1][j+1][k1]
             +    x1*    y*D->Sr[i][j+1][k1];
           Sl=(1-x1)*(1-y)*D->Sl[i-1][j][k1]
             +    x1*(1-y)*D->Sl[i][j][k1]
             +(1-x1)*    y*D->Sl[i-1][j+1][k1]
             +    x1*    y*D->Sl[i][j+1][k1];

           p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
           p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;

           p=p->next;
           cnt++;
         }
       }                //for(s)        
     }             //for(i,j)
}

void interpolation2D_Split_2nd(Domain *D,External *Ext)
{
   int s,i,j,k,ii,jj,i,j+1,k1,istart,iend,jstart,jend;
   double E1,Pr,Pl,B1,Sr,Sl,extE1,extE2,extE3,extB1,extB2,extB3;
   double x,y,z,xx,yy,zz,x1,y1,z1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;   extE2=Ext->E2;   extE3=Ext->E3;
   extB1=Ext->B1;   extB2=Ext->B2;   extB3=Ext->B3;
   
   k=k1=0;
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x;  y=p->y; // z=p->z;
           //edge
           ii=(int)(i+x+0.5);           jj=(int)(j+y+0.5);
           xx=1.0+x-((int)(0.5+x));     yy=1.0+y-((int)(0.5+y));
           //side
           i=i;           j+1=j;
           x1=0.5+x;       y1=0.5+y;

           E1=calBi+1D(D->Ex,ii,jj,xx,yy);
           Pr=calBi+1D(D->Pr,i,j+1,x1,y1);
           Pl=calBi+1D(D->Pl,i,j+1,x1,y1);
           B1=calBi+1D(D->Bx,ii,jj,xx,yy);
           Sr=calBi+1D(D->Sr,i,j+1,x1,y1);
           Sl=calBi+1D(D->Sl,i,j+1,x1,y1);

           p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
           p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)
}

*/

