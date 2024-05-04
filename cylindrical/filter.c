#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "mpi.h"

void MPI_filterNew_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share,int m);
void MPI_filterNew_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share,int m);
void MPI_filterCenter_Xminus(Domain *D,int numField,double ***f,int share);
void MPI_filterCenter_Xplus(Domain *D,int numField,double ***f,int share);


void filter_center(Domain *D,int numField,int iteration)
{
  int i,j,m,istart,iend,jstart,jend,numMode,iter;
  int nxSub,nySub,n,ii,jj; 
  double ***field,wz[3],sum;
  int myrank,nTasks;  

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;
  nxSub=D->nxSub+5;   nySub=D->nySub+5;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  field=(double ***)malloc(numField*sizeof(double **));
  for(n=0; n<numField; n++) {
    field[n]=(double **)malloc(numMode*sizeof(double *)); 
    for(m=0; m<numMode; m++) 
      field[n][m]=(double *)malloc(nxSub*sizeof(double *)); 
  }
      
  for(n=0; n<numField; n++)
    for(m=0; m<numMode; m++)
      for(i=0; i<nxSub; i++)
        field[n][m][i]=0.0;


  wz[0]=0.25; wz[1]=0.5; wz[2]=0.25;
  j=jstart;
  for(n=0; n<numField; n++)
    for(m=0; m<numMode; m++)
      for(i=istart; i<iend; i++) {
         sum=0.0;
         for(ii=0; ii<3; ii++) sum+=wz[ii]*D->shareF[n][m][i-1+ii][j];
         field[n][m][i]=sum;
      }

  if(D->L>1)  {
     MPI_filterCenter_Xplus(D,numField,field,3);
     MPI_filterCenter_Xminus(D,numField,field,3);
  } else ;

  m=0;
  for(n=0; n<numField; n++)
      for(i=0; i<nxSub; i++) 
         D->shareF[n][m][i][j]=field[n][m][i];
  for(n=0; n<numField; n++)
    for(m=2; m<numMode; m++)
      for(i=0; i<nxSub; i++) 
         D->shareF[n][m][i][j]=field[n][m][i];

  for(n=0; n<numField; n++) {
    for(m=0; m<numMode; m++) free(field[n][m]);
    free(field[n]);
  }
  free(field); 
}


/*
void filter_center(Domain *D,double ***dataR,double ***dataI,int ny)
{
  int i,j,m,istart,iend,jstart,jend,numMode,iter;
  int nxSub,nySub,n,ii,jj; 
  double ***valR,***valI;
  double sumR,sumI,wz[3];
  int myrank,nTasks;  

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;
  nxSub=D->nxSub+5;   nySub=D->nySub+5;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  valR=(double ***)malloc(numMode*sizeof(double ** ));
  valI=(double ***)malloc(numMode*sizeof(double ** ));
  for(m=0; m<numMode; m++) {
    valR[m]=(double **)malloc(nxSub*sizeof(double * ));
    valI[m]=(double **)malloc(nxSub*sizeof(double * ));
    for(i=0; i<nxSub; i++) {
      valR[m][i]=(double *)malloc(ny*sizeof(double  ));
      valI[m][i]=(double *)malloc(ny*sizeof(double  ));
	  }
  }
  for(m=0; m<numMode; m++)
    for(i=0; i<nxSub; i++)
      for(j=0; j<ny; j++) {
        valR[m][i][j]=0.0;
        valI[m][i][j]=0.0;
      }

  wz[0]=0.25; wz[1]=0.5; wz[2]=0.25;

  //first
  m=0;
  for(j=0; j<ny; j++) 
    for(i=istart; i<iend; i++) {
      sumR=0.0;
      for(ii=0; ii<3; ii++)
        sumR+=wz[ii]*dataR[m][i-1+ii][j+jstart];
      valR[m][i][j]=sumR;
    }
  for(m=1; m<numMode; m++) 
    for(j=0; j<ny; j++) 
      for(i=istart; i<iend; i++) {
        sumR=0.0;
        sumI=0.0;
        for(ii=0; ii<3; ii++) {
          sumR+=wz[ii]*dataR[m][i-1+ii][j+jstart];
          sumI+=wz[ii]*dataI[m][i-1+ii][j+jstart];
		    }
        valR[m][i][j]=sumR;
        valI[m][i][j]=sumI;
      }

  if(D->L>1)  {
    MPI_filterCenter_Xplus(D,valR,valI,ny,3);
    MPI_filterCenter_Xminus(D,valR,valI,ny,3);
  } else ;

  //second
  m=0;
  for(j=0; j<ny; j++) 
    for(i=istart; i<iend; i++) {
      sumR=0.0;
      for(ii=0; ii<3; ii++)
        sumR+=wz[ii]*valR[m][i-1+ii][j];
      dataR[m][i][j+jstart]=sumR;
    }		  

  for(m=1; m<numMode; m++) 
    for(j=0; j<ny; j++) 
      for(i=istart; i<iend; i++) {
        sumR=0.0;
        sumI=0.0;
        for(ii=0; ii<3; ii++) {
          sumR+=wz[ii]*valR[m][i-1+ii][j];
          sumI+=wz[ii]*valI[m][i-1+ii][j];
	      }
        dataR[m][i][j+jstart]=sumR;
        dataI[m][i][j+jstart]=sumI;
      }

  if(D->L>1)  {
    MPI_Transfer2F_Xplus(D,dataR,dataI,D->nySub+5,3);
    MPI_Transfer2F_Xminus(D,dataR,dataI,D->nySub+5,3);
  }  else ;

  for(m=0; m<numMode; m++) {
    for(i=0; i<nxSub; i++) { 
      free(valR[m][i]);
      free(valI[m][i]);
    }
    free(valR[m]); 
    free(valI[m]); 
  }
  free(valR); 
  free(valI); 
}
*/


void filter(Domain *D,double ***dataR,double ***dataI)
{
    int i,j,m,istart,iend,jstart,jend,numMode,iter;
    int nxSub,nySub,n,ii,jj; 
     double ***valR,***valI;
    double sumR,sumI,wz[3],wr[3];
    int myrank,nTasks;  

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;
    nxSub=D->nxSub+5;    nySub=D->nySub+5;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    valR=(double ***)malloc(1*sizeof(double ** ));
    valI=(double ***)malloc(1*sizeof(double ** ));
    for(m=0; m<1; m++) {
      valR[m]=(double **)malloc(nxSub*sizeof(double * ));
      valI[m]=(double **)malloc(nxSub*sizeof(double * ));
      for(i=0; i<nxSub; i++) {
        valR[m][i]=(double *)malloc(nySub*sizeof(double  ));
        valI[m][i]=(double *)malloc(nySub*sizeof(double  ));
		}
    }
    for(i=0; i<nxSub; i++)
      for(j=0; j<nySub; j++) {
         valR[0][i][j]=0.0;
         valI[0][i][j]=0.0;
      }

    wz[0]=0.25; wz[1]=0.5; wz[2]=0.25;

    m=0;
        //first
        j=jstart;
        wr[0]=0.5; wr[1]=0.5;
        for(i=istart; i<iend; i++) {
          sumR=0.0;
          for(ii=0; ii<3; ii++)
            for(jj=0; jj<2; jj++)
              sumR+=wz[ii]*wr[jj]*dataR[m][i-1+ii][j+jj];
          valR[0][i][j]=sumR;
        }

        wr[0]=0.25; wr[1]=0.5; wr[2]=0.25;
        for(i=istart; i<iend; i++)
          for(j=jstart+1; j<jend; j++) {
            sumR=0.0;
            for(ii=0; ii<3; ii++)
              for(jj=0; jj<3; jj++)
                sumR+=wz[ii]*wr[jj]*dataR[m][i-1+ii][j-1+jj];
            valR[0][i][j]=sumR;
          }

        if(D->L>1)  {
          MPI_filterNew_Xplus(D,valR,valI,nySub,3,0);
          MPI_filterNew_Xminus(D,valR,valI,nySub,3,0);
        } else ;

        //second
        j=jstart;
        wr[0]=0.5; wr[1]=0.5;
        for(i=istart; i<iend; i++) {
          sumR=0.0;
          for(ii=0; ii<3; ii++)
            for(jj=0; jj<2; jj++)
              sumR+=wz[ii]*wr[jj]*valR[0][i-1+ii][j+jj];
          dataR[m][i][j]=sumR;
        }

        wr[0]=0.25; wr[1]=0.5; wr[2]=0.25;
        for(i=istart; i<iend; i++)
          for(j=jstart+1; j<jend; j++) {
            sumR=0.0;
            for(ii=0; ii<3; ii++)
              for(jj=0; jj<3; jj++)
                sumR+=wz[ii]*wr[jj]*valR[0][i-1+ii][j-1+jj];
            dataR[m][i][j]=sumR;
          }

        if(D->L>1)  {
          MPI_filterNew_Xplus(D,dataR,dataI,nySub,3,m);
          MPI_filterNew_Xminus(D,dataR,dataI,nySub,3,m);
        } else ;
		  
    
    for(m=1; m<numMode; m++) {
        //first
        j=jstart;
        wr[0]=0.5; wr[1]=0.0;
        for(i=istart; i<iend; i++) {
          sumR=0.0;
          sumI=0.0;
          for(ii=0; ii<3; ii++)
            for(jj=0; jj<2; jj++) {
              sumR+=wz[ii]*wr[jj]*dataR[m][i-1+ii][j+jj];
              sumI+=wz[ii]*wr[jj]*dataI[m][i-1+ii][j+jj];
				}
          valR[0][i][j]=sumR;
          valI[0][i][j]=sumI;
        }

        wr[0]=0.25; wr[1]=0.5; wr[2]=0.25;
        for(i=istart; i<iend; i++)
          for(j=jstart+1; j<jend; j++) {
            sumR=0.0;
            sumI=0.0;
            for(ii=0; ii<3; ii++)
              for(jj=0; jj<3; jj++) {
                sumR+=wz[ii]*wr[jj]*dataR[m][i-1+ii][j-1+jj];
                sumI+=wz[ii]*wr[jj]*dataI[m][i-1+ii][j-1+jj];
				  }
            valR[0][i][j]=sumR;
            valI[0][i][j]=sumI;
          }

        if(D->L>1)  {
          MPI_filterNew_Xplus(D,valR,valI,nySub,3,0);
          MPI_filterNew_Xminus(D,valR,valI,nySub,3,0);
        } else ;

        //second
        j=jstart;
        wr[0]=0.5; wr[1]=0.0;
        for(i=istart; i<iend; i++) {
          sumR=0.0;
          sumI=0.0;
          for(ii=0; ii<3; ii++)
            for(jj=0; jj<2; jj++) {
              sumR+=wz[ii]*wr[jj]*valR[0][i-1+ii][j+jj];
              sumI+=wz[ii]*wr[jj]*valI[0][i-1+ii][j+jj];
				}
          dataR[m][i][j]=sumR;
          dataI[m][i][j]=sumI;
        }

        wr[0]=0.25; wr[1]=0.5; wr[2]=0.25;
        for(i=istart; i<iend; i++)
          for(j=jstart+1; j<jend; j++) {
            sumR=0.0;
            sumI=0.0;
            for(ii=0; ii<3; ii++)
              for(jj=0; jj<3; jj++) {
                sumR+=wz[ii]*wr[jj]*valR[0][i-1+ii][j-1+jj];
                sumI+=wz[ii]*wr[jj]*valI[0][i-1+ii][j-1+jj];
				  }
            dataR[m][i][j]=sumR;
            dataI[m][i][j]=sumI;
          }

        if(D->L>1)  {
          MPI_filterNew_Xplus(D,dataR,dataI,nySub,3,m);
          MPI_filterNew_Xminus(D,dataR,dataI,nySub,3,m);
        } else ;

    }
	
    for(i=0; i<nxSub; i++) { 
      free(valR[0][i]);
      free(valI[0][i]);
    }
    free(valR[0]); free(valR); 
    free(valI[0]); free(valI); 
}

void MPI_filterCenter_Xminus(Domain *D,int numField,double ***f,int share)
{
    int i,j,n,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank;
    double *data;

    MPI_Status status;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    rank=myrank/D->M;

    num=numField*numMode*share;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);
      start=0;
      for(n=0; n<numField; n++) 
        for(m=0; m<numMode; m++) {
          for(i=0; i<share; i++) f[n][m][iend+i]=data[start+i]; 
          start+=share;
        }
    }
    else if(rank%2==1) {
      start=0;
      for(n=0; n<numField; n++) 
        for(m=0; m<numMode; m++) {
          for(i=0; i<share; i++) data[start+i]=f[n][m][istart+i];  
          start+=share;
	     }
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);
      start=0;
      for(n=0; n<numField; n++) 
        for(m=0; m<numMode; m++) {
          for(i=0; i<share; i++) f[n][m][iend+i]=data[start+i]; 
          start+=share;
        }
    }
    else if(rank%2==0 && rank!=0) {
      start=0;
      for(n=0; n<numField; n++) 
        for(m=0; m<numMode; m++) {
          for(i=0; i<share; i++) data[start+i]=f[n][m][istart+i];  
          start+=share;
	     }
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_filterCenter_Xplus(Domain *D,int numField,double ***f,int share)
{
    int i,n,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start;
    double *data;

    MPI_Status status;

    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    rank=myrank/D->M;

    num=numField*numMode*(share-1);
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    if(rank%2==1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);
      start=0;
      for(n=0; n<numField; n++) 
        for(m=0; m<numMode; m++) {   
          for(i=1; i<share; i++) f[n][m][istart-i]=data[start+(i-1)];
          start+=(share-1);
        }
    }
    else if(rank%2==0 && rank!=D->L-1) {
      start=0;
      for(n=0; n<numField; n++) 
        for(m=0; m<numMode; m++)  { 
          for(i=1; i<share; i++)  data[start+(i-1)]=f[n][m][iend-i];
          start+=(share-1);
        }
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             

    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
      start=0;
      for(n=0; n<numField; n++) 
        for(m=0; m<numMode; m++) {   
          for(i=1; i<share; i++) f[n][m][istart-i]=data[start+i-1];
          start+=(share-1);
        }
    }
    else if(rank%2==1 && rank!=D->L-1) {
      start=0;
      for(n=0; n<numField; n++) 
        for(m=0; m<numMode; m++)  { 
          for(i=1; i<share; i++)  data[start+i-1]=f[n][m][iend-i];
          start+=(share-1);
        }
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
	 
    free(data);
}


void MPI_filterNew_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share,int m)
{
    int i,j,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank;
    double *data;

    MPI_Status status;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    rank=myrank/D->M;

    num=2*ny*3;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
    for(i=0; i<share; i++) {
      for(j=0; j<ny; j++)  data[start+j]=f1[m][i+istart][j];  start+=ny;
      for(j=0; j<ny; j++)  data[start+j]=f2[m][i+istart][j];  start+=ny;
	 }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);
      start=0;
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
      }

    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);
      start=0;
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_filterNew_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share,int m)
{
    int i,j,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start;
    double *data;

    MPI_Status status;

    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    rank=myrank/D->M;

    num=2*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
      }

    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);
      start=0;
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
      }

    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
      start=0;
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}
