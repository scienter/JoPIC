#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"
 

void restoreDumpHDF(Domain *D,int iteration)
{
    int i,j,k,s,n,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int nxSub,nySub,nzSub,nSpecies,totalCnt;
    int *offSet,*counts;
    int *dataI,*dataJ,*dataK,*dataCore,*dataIndex;
    double *dataX,*dataY,*dataZ,*dataPx,*dataPy,*dataPz;
    int offset[3];
    char name[100],name2[100];
    void restoreIntMeta();
    void restoreField2D();
    void restoreField3D();
    void restoreIntData();
    void restoreDoubleData();
    void restoreIntArray();
    ptclList *p;
    Particle ***particle;
    particle=D->particle;

    int myrank, nTasks;    
    MPI_Status mpi_status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,filespace;
    hsize_t metaDim[1];
    herr_t status;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;
    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    metaDim[0]=1;

    offSet=(int *)malloc(nTasks*sizeof(int ));
    counts=(int *)malloc(nTasks*sizeof(int ));

    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Split-1)*3+2:
      ny=D->ny+5;      
      nx=D->nx+5;     
      istart=0;
      iend+=3; 
      jstart=0;
      jend+=3;
      nxSub+=5;
      nySub+=5;
      offset[0]=0;
      offset[1]=D->minYSub-D->minYDomain;
      
      //save metaData for field
      sprintf(name,"dump%d.h5",iteration);
      if(myrank==0)
        restoreIntMeta(name,"/minXSub",&(D->minXSub));
      else	;
      MPI_Bcast(&(D->minXSub),1,MPI_INT,0,MPI_COMM_WORLD);

      restoreField2D(D->Ex,name,"/Ex",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->Pr,name,"/Pr",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->Pl,name,"/Pl",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->Bx,name,"/Bx",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->Sr,name,"/Sr",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->Sl,name,"/Sl",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->ExC,name,"/ExC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->PrC,name,"/PrC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->PlC,name,"/PlC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->BxC,name,"/BxC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->SrC,name,"/SrC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->SlC,name,"/SlC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->Jx,name,"/Jx",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->Jy,name,"/Jy",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->Jz,name,"/Jz",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->JxOld,name,"/JxOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->JyOld,name,"/JyOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      restoreField2D(D->JzOld,name,"/JzOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);

      //restore metaData for particle
      if(myrank==0)
        restoreIntMeta(name,"/nSpecies",&nSpecies);
      else	;
      MPI_Bcast(&nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);

      //restore particle
      for(s=0; s<nSpecies; s++)
      {
        sprintf(name2,"%dParticle/totalCnt",s);
        if(myrank==0)
          restoreIntMeta(name,name2,&totalCnt);
        else	;
        MPI_Bcast(&totalCnt,1,MPI_INT,0,MPI_COMM_WORLD);
        sprintf(name2,"%dParticle/offSet",s);
        if(myrank==0)
          restoreIntArray(name,name2,offSet,nTasks);
        else	;
        MPI_Bcast(offSet,nTasks,MPI_INT,0,MPI_COMM_WORLD);
        sprintf(name2,"%dParticle/counts",s);
        if(myrank==0)
          restoreIntArray(name,name2,counts,nTasks);
        else	;
        MPI_Bcast(counts,nTasks,MPI_INT,0,MPI_COMM_WORLD);

        if(counts[myrank]>0)
        {
          dataI = (int *)malloc(counts[myrank]*sizeof(int ));
          dataJ = (int *)malloc(counts[myrank]*sizeof(int ));
          dataCore = (int *)malloc(counts[myrank]*sizeof(int ));
          dataIndex = (int *)malloc(counts[myrank]*sizeof(int ));
          dataX = (double *)malloc(counts[myrank]*sizeof(double ));
          dataY = (double *)malloc(counts[myrank]*sizeof(double ));
          dataPx = (double *)malloc(counts[myrank]*sizeof(double ));
          dataPy = (double *)malloc(counts[myrank]*sizeof(double ));
          dataPz = (double *)malloc(counts[myrank]*sizeof(double ));
          sprintf(name2,"%dParticle/i",s);
          restoreIntData(name,name2,dataI,totalCnt,counts[myrank],offSet[myrank]);
          MPI_Barrier(MPI_COMM_WORLD);
          sprintf(name2,"%dParticle/j",s);
          restoreIntData(name,name2,dataJ,totalCnt,counts[myrank],offSet[myrank]);
          MPI_Barrier(MPI_COMM_WORLD);
          sprintf(name2,"%dParticle/core",s);
          restoreIntData(name,name2,dataCore,totalCnt,counts[myrank],offSet[myrank]);
          MPI_Barrier(MPI_COMM_WORLD);
          sprintf(name2,"%dParticle/index",s);
          restoreIntData(name,name2,dataIndex,totalCnt,counts[myrank],offSet[myrank]);
          MPI_Barrier(MPI_COMM_WORLD);
          sprintf(name2,"%dParticle/x",s);
          restoreDoubleData(name,name2,dataX,totalCnt,counts[myrank],offSet[myrank]);
          MPI_Barrier(MPI_COMM_WORLD);
          sprintf(name2,"%dParticle/y",s);
          restoreDoubleData(name,name2,dataY,totalCnt,counts[myrank],offSet[myrank]);
          MPI_Barrier(MPI_COMM_WORLD);
          sprintf(name2,"%dParticle/px",s);
          restoreDoubleData(name,name2,dataPx,totalCnt,counts[myrank],offSet[myrank]);
          MPI_Barrier(MPI_COMM_WORLD);
          sprintf(name2,"%dParticle/py",s);
          restoreDoubleData(name,name2,dataPy,totalCnt,counts[myrank],offSet[myrank]);
          MPI_Barrier(MPI_COMM_WORLD);
          sprintf(name2,"%dParticle/pz",s);
          restoreDoubleData(name,name2,dataPz,totalCnt,counts[myrank],offSet[myrank]);
          MPI_Barrier(MPI_COMM_WORLD);

          for(n=0; n<counts[myrank]; n++)
          {
            i=dataI[n];
            j=dataJ[n];
            p = (ptclList *)malloc(sizeof(ptclList));
            p->next = particle[i][j][0].head[s]->pt;
            particle[i][j][0].head[s]->pt=p;
            
            p->x=dataX[n];
            p->y=dataY[n];
            p->z=0.0;
            p->p1=dataPx[n];
            p->p2=dataPy[n];
            p->p3=dataPz[n];
            p->index=dataIndex[n];
            p->core=dataCore[n];
          }
          free(dataI);
          free(dataJ);
          free(dataCore);
          free(dataIndex);
          free(dataX);
          free(dataY);
          free(dataPx);
          free(dataPy);
          free(dataPz);
        }
        else	;
      }			//End of for(nSpecies)
    break;

    //3D
    case (Split-1)*3+3:
      ny=D->ny+5;      
      nx=D->nx+5;     
      nz=D->nz+5;      
      istart=0;
      iend+=3; 
      jstart=0;
      jend+=3;
      kstart=0;
      kend+=3;
      nxSub+=5;
      nySub+=5;
      nzSub+=5;
      offset[1]=0;
      offset[0]=D->minYSub-D->minYDomain;
      offset[2]=D->minZSub-D->minZDomain;
      
      //save metaData for field
      sprintf(name,"dump%d.h5",iteration);
      if(myrank==0)
        restoreIntMeta(name,"/minXSub",&(D->minXSub));
      else	;
      MPI_Bcast(&(D->minXSub),1,MPI_INT,0,MPI_COMM_WORLD);
//here
      restoreField3D(D->Ex,name,"/Ex",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Pr,name,"/Pr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Pl,name,"/Pl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Bx,name,"/Bx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Sr,name,"/Sr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Sl,name,"/Sl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->ExC,name,"/ExC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->PrC,name,"/PrC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->PlC,name,"/PlC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->BxC,name,"/BxC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->SrC,name,"/SrC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->SlC,name,"/SlC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Jx,name,"/Jx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Jy,name,"/Jy",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Jz,name,"/Jz",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->JxOld,name,"/JxOld",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->JyOld,name,"/JyOld",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->JzOld,name,"/JzOld",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);

      //restore metaData for particle
      if(myrank==0)
        restoreIntMeta(name,"/nSpecies",&nSpecies);
      else	;
      MPI_Bcast(&nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);

      //restore particle
      for(s=0; s<nSpecies; s++)
      {
        sprintf(name2,"%dParticle/totalCnt",s);
        if(myrank==0)
          restoreIntMeta(name,name2,&totalCnt);
        else	;
        MPI_Bcast(&totalCnt,1,MPI_INT,0,MPI_COMM_WORLD);
        sprintf(name2,"%dParticle/offSet",s);
        if(myrank==0)
          restoreIntArray(name,name2,offSet,nTasks);
        else	;
        MPI_Bcast(offSet,nTasks,MPI_INT,0,MPI_COMM_WORLD);
        sprintf(name2,"%dParticle/counts",s);
        if(myrank==0)
          restoreIntArray(name,name2,counts,nTasks);
        else	;
        MPI_Bcast(counts,nTasks,MPI_INT,0,MPI_COMM_WORLD);

        if(counts[myrank]>0)
        {
          dataI = (int *)malloc(counts[myrank]*sizeof(int ));
          dataJ = (int *)malloc(counts[myrank]*sizeof(int ));
          dataK = (int *)malloc(counts[myrank]*sizeof(int ));
          dataCore = (int *)malloc(counts[myrank]*sizeof(int ));
          dataIndex = (int *)malloc(counts[myrank]*sizeof(int ));
          dataX = (double *)malloc(counts[myrank]*sizeof(double ));
          dataY = (double *)malloc(counts[myrank]*sizeof(double ));
          dataZ = (double *)malloc(counts[myrank]*sizeof(double ));
          dataPx = (double *)malloc(counts[myrank]*sizeof(double ));
          dataPy = (double *)malloc(counts[myrank]*sizeof(double ));
          dataPz = (double *)malloc(counts[myrank]*sizeof(double ));
          sprintf(name2,"%dParticle/i",s);
          restoreIntData(name,name2,dataI,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/j",s);
          restoreIntData(name,name2,dataJ,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/k",s);
          restoreIntData(name,name2,dataK,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/core",s);
          restoreIntData(name,name2,dataCore,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/index",s);
          restoreIntData(name,name2,dataIndex,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/x",s);
          restoreDoubleData(name,name2,dataX,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/y",s);
          restoreDoubleData(name,name2,dataY,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/z",s);
          restoreDoubleData(name,name2,dataZ,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/px",s);
          restoreDoubleData(name,name2,dataPx,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/py",s);
          restoreDoubleData(name,name2,dataPy,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/pz",s);
          restoreDoubleData(name,name2,dataPz,totalCnt,counts[myrank],offSet[myrank]);

          for(n=0; n<counts[myrank]; n++)
          {
            i=dataI[n];
            j=dataJ[n];
            k=dataK[n];
            p = (ptclList *)malloc(sizeof(ptclList));
            p->next = particle[i][j][k].head[s]->pt;
            particle[i][j][k].head[s]->pt=p;
            
            p->x=dataX[n];
            p->y=dataY[n];
            p->z=dataZ[n];
            p->p1=dataPx[n];
            p->p2=dataPy[n];
            p->p3=dataPz[n];
            p->index=dataIndex[n];
            p->core=dataCore[n];
          }

          free(dataI);
          free(dataJ);
          free(dataK);
          free(dataCore);
          free(dataIndex);
          free(dataX);
          free(dataY);
          free(dataZ);
          free(dataPx);
          free(dataPy);
          free(dataPz);
        }
      }			//End of for(nSpecies)
    break;

    free(offSet);
    free(counts);
  }		//End of switch(dimension....)
}

void restoreIntData(char *fileName,char *dataName,int *data,int totalCnt,int dataNum,int offSet)
{
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id,group_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[1],count[1],offset[1];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=totalCnt;     
  filespace=H5Screate_simple(1,dimsf,NULL);

  count[0]=dataNum;
  offset[0]=offSet;
  memspace=H5Screate_simple(1,count,NULL);

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
  status = H5Dread(dset_id, H5T_NATIVE_INT,memspace,subfilespace,plist_id,data);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restoreDoubleData(char *fileName,char *dataName,double *data,int totalCnt,int dataNum,int offSet)
{
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id,group_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[1],count[1],offset[1];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=totalCnt;     
  filespace=H5Screate_simple(1,dimsf,NULL);

  count[0]=dataNum;
  offset[0]=offSet;
  memspace=H5Screate_simple(1,count,NULL);

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,data);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restoreField2D(double ***data,char *fileName,char *dataName,int nx,int ny,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int *offSet)
{
  int i,j,k,start;
  double *field;
  char name[100];
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[3],count[3],offset[3];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//  H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//  MPI_Barrier(MPI_COMM_WORLD);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=nx;      
  dimsf[1]=ny;     
  filespace=H5Screate_simple(2,dimsf,NULL);

  count[1]=nySub;
  count[0]=nxSub;
  offset[1]=offSet[1];
  offset[0]=offSet[0];
  memspace=H5Screate_simple(2,count,NULL);

  field = (double *)malloc(nxSub*nySub*sizeof(double ));

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
  start=0;
  for(i=istart; i<iend; i++)
  {
    for(j=jstart; j<jend; j++)
      data[i][j][0]=field[start+j-jstart];
    start+=nySub;
  }
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
  free(field);
}

void restoreField3D(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet)
{
  int i,j,k,start;
  double *field;
  char name[100];
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[3],count[3],offset[3];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//  H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//  MPI_Barrier(MPI_COMM_WORLD);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=ny;      
  dimsf[1]=nx;     
  dimsf[2]=nz;     
  filespace=H5Screate_simple(3,dimsf,NULL);

  count[0]=nySub;
  count[1]=nxSub;
  count[2]=nzSub;
  offset[1]=offSet[1];
  offset[0]=offSet[0];
  offset[2]=offSet[2];
  memspace=H5Screate_simple(3,count,NULL);

  field = (double *)malloc(nxSub*nySub*nzSub*sizeof(double ));

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
  start=0;
  for(j=jstart; j<jend; j++)
    for(i=istart; i<iend; i++)
    {
      for(k=kstart; k<kend; k++)
        data[i][j][k]=field[start+k-kstart];
      start+=nzSub;
    }
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
  free(field);
}

void restoreIntMeta(char *fileName,char *dataName,int *data)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restoreIntArray(char *fileName,char *dataName,int *data,int nTasks)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[2];
  herr_t status;

  metaDim[0]=nTasks;
  metaDim[1]=1;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(2,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

