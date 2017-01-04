#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"
#include "hdf5.h"
#include "constants.h"
#include "math.h"

void saveTracking(Domain *D)
{
  int n,i,dataNum,numberData;
  double *transfer,*hdfField;
  int myrank,nTasks,rnk,totalCnt;
  double xx,yy,zz,px,py,pz,id,core,step;
  char name[100],dataName[100],fileName[100];
  FILE *out;

  hid_t file_id;
  void saveDoubleArray(char *fileName,char *dataName,double *data,int totalCnt);

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  numberData=(D->maxStep-D->trackStart)/D->trackSaveStep+1;
  dataNum=9*numberData;
  transfer=(double *)malloc(dataNum*sizeof(double));
  hdfField=(double *)malloc(numberData*sizeof(double));


  for(n=0; n<D->idNums; n++)
  {
    for(i=0; i<dataNum; i++)
      transfer[i]=0.0;

    for(i=0; i<numberData; i++)
    {
      transfer[i*9+0]=D->track[n][i].x;
      transfer[i*9+1]=D->track[n][i].y;
      transfer[i*9+2]=D->track[n][i].z;
      transfer[i*9+3]=D->track[n][i].px;
      transfer[i*9+4]=D->track[n][i].py;
      transfer[i*9+5]=D->track[n][i].pz;
      transfer[i*9+6]=(double)(D->track[n][i].id);
      transfer[i*9+7]=(double)(D->track[n][i].core);
      transfer[i*9+8]=(double)(D->track[n][i].step);
    }

    if(myrank==0)
    {
      for(rnk=1; rnk<nTasks; rnk++)
      {
        MPI_Recv(transfer,dataNum,MPI_DOUBLE,rnk,rnk,MPI_COMM_WORLD,&status);
        for(i=0; i<numberData; i++)
        {
          D->track[n][i].x+=transfer[i*9+0];
          D->track[n][i].y+=transfer[i*9+1];
          D->track[n][i].z+=transfer[i*9+2];
          D->track[n][i].px+=transfer[i*9+3];
          D->track[n][i].py+=transfer[i*9+4];
          D->track[n][i].pz+=transfer[i*9+5];
          D->track[n][i].id+=(int)(transfer[i*9+6]);
          D->track[n][i].core+=(int)(transfer[i*9+7]);
          D->track[n][i].step+=(int)(transfer[i*9+8]);
        }
      }
    }
    else
      MPI_Send(transfer,dataNum,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if(myrank==0)
    {
      if(D->saveMode==TXT)
      {
        sprintf(name,"%dTrack%d_%d",D->trackS[n],D->trackID[n],D->trackCore[n]);
        out=fopen(name,"w");
        for(i=0; i<numberData; i++)
        {
          xx=D->track[n][i].x;
          yy=D->track[n][i].y;
          zz=D->track[n][i].z;
          px=D->track[n][i].px;
          py=D->track[n][i].py;
          pz=D->track[n][i].pz;
          id=D->track[n][i].id;
          core=D->track[n][i].core;
          step=D->track[n][i].step;
          if(xx>0)
            fprintf(out,"%.10g %.10g %.10g %.10g %.10g %.10g %g %g %g\n",xx,yy,zz,px,py,pz,id,core,step);
          else	;
        }
        fclose(out);
      }
      else if(D->saveMode==HDF)
      {
        sprintf(fileName,"%dTrack%d_%d.h5",D->trackS[n],D->trackID[n],D->trackCore[n]);
        file_id=H5Fcreate(fileName,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        H5Fclose(file_id);
        totalCnt=numberData;

        sprintf(dataName,"x");
        for(i=0; i<numberData; i++)
          hdfField[i]=D->track[n][i].x;
        saveDoubleArray(fileName,dataName,hdfField,totalCnt);
        sprintf(dataName,"y");
        for(i=0; i<numberData; i++)
          hdfField[i]=D->track[n][i].y;
        saveDoubleArray(fileName,dataName,hdfField,totalCnt);
        sprintf(dataName,"z");
        for(i=0; i<numberData; i++)
          hdfField[i]=D->track[n][i].z;
        saveDoubleArray(fileName,dataName,hdfField,totalCnt);
        sprintf(dataName,"px");
        for(i=0; i<numberData; i++)
          hdfField[i]=D->track[n][i].px;
        saveDoubleArray(fileName,dataName,hdfField,totalCnt);
        sprintf(dataName,"py");
        for(i=0; i<numberData; i++)
          hdfField[i]=D->track[n][i].py;
        saveDoubleArray(fileName,dataName,hdfField,totalCnt);
        sprintf(dataName,"pz");
        for(i=0; i<numberData; i++)
          hdfField[i]=D->track[n][i].pz;
        saveDoubleArray(fileName,dataName,hdfField,totalCnt);
        sprintf(dataName,"step");
        for(i=0; i<numberData; i++)
          hdfField[i]=D->track[n][i].step;
        saveDoubleArray(fileName,dataName,hdfField,totalCnt);
        sprintf(dataName,"totalCnt");
        saveIntMeta(fileName,dataName,&totalCnt);

        printf("%s is made.\n",fileName);
      }
      else	;
    }
    else 	; //End of (myrank==0)

  }	//for(n)
  free(transfer);
  free(hdfField);

}


void trackID(Domain *D,int iteration,int istart,int iend,int jstart,int jend,int kstart,int kend)
{
  int i,j,k,n,iter,s,ii,index,transfer;
  double wp,kp;
  Particle ***particle;
  particle=D->particle;
  ptclList *p;

//  wp=sqrt(D->density*eCharge*eCharge/eMass/eps0);
//  kp=wp/velocityC;

  index=(iteration-D->trackStart)/D->trackSaveStep;

  for(n=0; n<D->idNums; n++)
  {
    transfer=0;
    s=D->trackS[n];
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {
            if(p->index==D->trackID[n] && p->core==D->trackCore[n])
            {
              D->track[n][index].x=(i-D->istart+D->minXSub+p->x)*D->dx*D->lambda;
              D->track[n][index].y=(j-D->jstart+D->minYSub+p->y)*D->dy*D->lambda;
              D->track[n][index].z=(k-D->kstart+D->minZSub+p->z)*D->dz*D->lambda;
              D->track[n][index].px=p->p1;
              D->track[n][index].py=p->p2;
              D->track[n][index].pz=p->p3;
              D->track[n][index].id=p->index;
              D->track[n][index].core=p->core;
              D->track[n][index].step=iteration;
              transfer=1;
            }
            p=p->next;
          }
        }
  }		//End of for(idNums)

}

void saveDoubleArray(char *fileName,char *dataName,double *data,int totalCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=totalCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

