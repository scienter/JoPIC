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
    int i,ii,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz,bias,offSetY;
    int nxSub,nySub,nzSub,nSpecies,totalCnt,tmp,cnt,remain,min,max,numData;
    int rank,index,indexI,indexJ,indexK,rnk,subrnk;
    int *cntOffSet,*dataI,*dataJ,*dataK,*dataCore,*dataIndex,*numInRank,*indexCnt;
    float *dataX,*dataY,*dataPx,*dataPy,*dataPz,**floatData,**recvData,*sendData;
    int *recvCnt,sendCnt,tmpInt;
    int offset[3];
    char name[100],name2[100];
    void restoreIntMeta();
    void restoreField2D();
    void restoreField3D();
    void restoreIntData();
    void restoreFloatData();
    ptclList *p;
    Particle ***particle;
    particle=D->particle;

    int myrank, nTasks;    
    MPI_Status mpi_status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    int minYSub[D->M],maxYSub[D->M],minZSub[D->N],maxZSub[D->N];

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

    cntOffSet=(int *)malloc(nTasks*sizeof(int ));
    recvCnt = (int *)malloc(nTasks*sizeof(int ));
    for(i=0; i<nTasks; i++)
      recvCnt[i]=0;

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

      //save minYSub,maxYSub,minZSub,maxZSub
      MPI_Gather(&(D->minYSub),1,MPI_INT,minYSub,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(minYSub,D->M,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Gather(&(D->maxYSub),1,MPI_INT,maxYSub,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(maxYSub,D->M,MPI_INT,0,MPI_COMM_WORLD);

      //restore metaData for particle
      if(myrank==0)
        restoreIntMeta(name,"/nSpecies",&nSpecies);
      else	;
      MPI_Bcast(&nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);

      restoreField2D(D->Ex,name,"/Ex",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->Pr,name,"/Pr",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->Pl,name,"/Pl",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->Bx,name,"/Bx",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->Sr,name,"/Sr",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->Sl,name,"/Sl",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->ExC,name,"/ExC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->PrC,name,"/PrC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->PlC,name,"/PlC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->BxC,name,"/BxC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->SrC,name,"/SrC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->SlC,name,"/SlC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->Jx,name,"/Jx",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->Jy,name,"/Jy",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->Jz,name,"/Jz",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->JxOld,name,"/JxOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->JyOld,name,"/JyOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      restoreField2D(D->JzOld,name,"/JzOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);

      //restore particle
      for(s=0; s<nSpecies; s++)
      {
        sprintf(name2,"%dParticle/totalCnt",s);
        if(myrank==0)
          restoreIntMeta(name,name2,&totalCnt);
        else	;
        MPI_Bcast(&totalCnt,1,MPI_INT,0,MPI_COMM_WORLD);
      
        cnt=totalCnt/nTasks;
        remain=totalCnt%nTasks;
        min=max=0;
        for(i=0; i<nTasks; i++)
        {
          if(i<remain) tmp=cnt+1;
          else	     tmp=cnt;
          min=max;
          max=min+tmp;
          if(i==myrank)
          {
            cntOffSet[i]=min;
            numData=tmp;
          }
        }

        if(numData>0)
        {
          dataI = (int *)malloc(numData*sizeof(int ));
          dataJ = (int *)malloc(numData*sizeof(int ));
          dataCore = (int *)malloc(numData*sizeof(int ));
          dataIndex = (int *)malloc(numData*sizeof(int ));
          dataX = (float *)malloc(numData*sizeof(float ));
          dataY = (float *)malloc(numData*sizeof(float ));
          dataPx = (float *)malloc(numData*sizeof(float ));
          dataPy = (float *)malloc(numData*sizeof(float ));
          dataPz = (float *)malloc(numData*sizeof(float ));
          numInRank=(int *)malloc(D->M*sizeof(int ));
          sprintf(name2,"%dParticle/i",s);

          restoreIntData(name,name2,dataI,totalCnt,numData,cntOffSet[myrank]);
          sprintf(name2,"%dParticle/j",s);
          restoreIntData(name,name2,dataJ,totalCnt,numData,cntOffSet[myrank]);
          sprintf(name2,"%dParticle/core",s);
          restoreIntData(name,name2,dataCore,totalCnt,numData,cntOffSet[myrank]);
          sprintf(name2,"%dParticle/index",s);
          restoreIntData(name,name2,dataIndex,totalCnt,numData,cntOffSet[myrank]);
          sprintf(name2,"%dParticle/x",s);
          restoreFloatData(name,name2,dataX,totalCnt,numData,cntOffSet[myrank]);
          sprintf(name2,"%dParticle/y",s);
          restoreFloatData(name,name2,dataY,totalCnt,numData,cntOffSet[myrank]);
          sprintf(name2,"%dParticle/px",s);
          restoreFloatData(name,name2,dataPx,totalCnt,numData,cntOffSet[myrank]);
          sprintf(name2,"%dParticle/py",s);
          restoreFloatData(name,name2,dataPy,totalCnt,numData,cntOffSet[myrank]);
          sprintf(name2,"%dParticle/pz",s);
          restoreFloatData(name,name2,dataPz,totalCnt,numData,cntOffSet[myrank]);

          for(i=0; i<D->M; i++)
            numInRank[i]=0;
          cnt=D->ny/D->M;
          remain=D->ny%D->M;
          for(i=0; i<numData; i++)
          {
            index=dataJ[i]-D->minYDomain;
            if(index<remain*(cnt+1))   rank=index/(cnt+1);
            else     rank=(index-remain*(cnt+1))/cnt+remain;
            numInRank[rank]+=1;
          }
          floatData=(float **)malloc(nTasks*sizeof(float *));
          recvData=(float **)malloc(nTasks*sizeof(float *));
          indexCnt=(int *)malloc(nTasks*sizeof(int ));
          for(rnk=0; rnk<nTasks; rnk++)
          {
            floatData[rnk]=(float *)malloc(numInRank[rnk]*8*sizeof(float ));
            indexCnt[rnk]=0;
          }

          for(i=0; i<numData; i++)
          {
            index=dataJ[i]-D->minYDomain;
            indexI=dataI[i]+D->istart;
            indexJ=dataJ[i]-D->minYSub+D->jstart;
            if(index<remain*(cnt+1))   rank=index/(cnt+1);
            else     rank=(index-remain*(cnt+1))/cnt+remain;
            if(rank==myrank)
            {
              p = (ptclList *)malloc(sizeof(ptclList));
              p->next = particle[indexI][indexJ][0].head[s]->pt;
              particle[indexI][indexJ][0].head[s]->pt=p;
            
              p->x=dataX[i];
              p->y=dataY[i];
              p->p1=dataPx[i];
              p->p2=dataPy[i];
              p->p3=dataPz[i];
              p->index=dataIndex[i];
              p->core=dataCore[i];
            }
            else
            {
              ii=indexCnt[rank]*8;
              floatData[rank][ii]=dataX[i]+indexI;
              floatData[rank][ii+1]=(float)(dataJ[i]);
              floatData[rank][ii+2]=dataY[i];
              floatData[rank][ii+3]=dataPx[i];
              floatData[rank][ii+4]=dataPy[i];
              floatData[rank][ii+5]=dataPz[i];
              floatData[rank][ii+6]=(float)(dataIndex[i]);
              floatData[rank][ii+7]=(float)(dataCore[i]);
              indexCnt[rank]+=1;
            }
          }	//End of for(dataNum)
          free(dataI);
          free(dataJ);
          free(dataCore);
          free(dataIndex);
          free(dataX);
          free(dataY);
          free(dataPx);
          free(dataPy);
          free(dataPz);

          //send and receive data count
          for(rnk=0; rnk<nTasks; rnk++)
            recvCnt[rnk]=0;

          for(rnk=0; rnk<nTasks; rnk++)
          {
            if(myrank!=rnk)
              MPI_Send(&indexCnt[rnk],1,MPI_INT,rnk,myrank,MPI_COMM_WORLD);
          }
          for(rnk=0; rnk<nTasks; rnk++)
          {
            if(myrank!=rnk)
              MPI_Recv(&recvCnt[rnk],1,MPI_INT,rnk,rnk,MPI_COMM_WORLD,&mpi_status);
          }          
          MPI_Barrier(MPI_COMM_WORLD);

          //send and receive data
          for(rnk=0; rnk<nTasks; rnk++)
          {
            if(myrank!=rnk)
            {
              tmpInt=indexCnt[rnk]*8;
              MPI_Send(floatData[rnk],tmpInt,MPI_FLOAT,rnk,myrank,MPI_COMM_WORLD);
              free(floatData[rnk]);
            }
            recvData[rnk]=(float *)malloc(recvCnt[rnk]*8*sizeof(float ));
          }          
          for(rnk=0; rnk<nTasks; rnk++)
          {
            if(myrank==rnk)
            {
              for(subrnk=0; subrnk<nTasks; subrnk++)
              {
                if(myrank!=subrnk)
                {
                  tmpInt=recvCnt[subrnk]*8;
                  MPI_Recv(recvData[subrnk],tmpInt,MPI_FLOAT,subrnk,subrnk,MPI_COMM_WORLD,&mpi_status);
                  for(i=0; i<recvCnt[subrnk]; i++)
                  {
                    ii=i*8;
                    indexI=(int)(recvData[subrnk][ii]);
                    indexJ=((int)(recvData[subrnk][ii+1]))-D->minYSub+D->jstart;
                    p = (ptclList *)malloc(sizeof(ptclList));
                    p->next = particle[indexI][indexJ][0].head[s]->pt;
                    particle[indexI][indexJ][0].head[s]->pt=p;
                    p->x=recvData[subrnk][ii]-indexI;
                    p->y=recvData[subrnk][ii+2];
                    p->p1=recvData[subrnk][ii+3];
                    p->p2=recvData[subrnk][ii+4];
                    p->p3=recvData[subrnk][ii+5];
                    p->index=(int)(recvData[subrnk][ii+6]);
                    p->core=(int)(recvData[subrnk][ii+7]);
                  }
                  free(recvData[subrnk]);
                }
              }		//End of for(subrnk)
            }
          }		//End of for(rnk)          
          free(numInRank);
        }	
        else 	; 	//End of if(dataNum>0)

      }			//End of for(nSpecies)

      

/*      FILE *out;
      float x,y,Pr;
      int i,j;
      sprintf(name,"field%d",myrank);
      out=fopen(name,"w");
      for(i=2; i<iend-3; i++)
      {
        for(j=2; j<jend-3; j++)
        {
          x=(i-2+D->minXSub)*D->dx*D->lambda;
          y=(j-2+D->minYSub)*D->dy*D->lambda;
          Pr=D->Pr[i][j][0];
          fprintf(out,"%g %g %g\n",x,y,Pr);
        }
        fprintf(out,"\n");
      }
*/
    break;

/*
    //3D
    case (Split-1)*3+3:
      dimsf[0]=D->ny;      
      dimsf[1]=D->nx;     
      dimsf[2]=D->nz;     
      nx=D->nx;
      ny=D->ny;
      nz=D->nz;
      filespace=H5Screate_simple(D->dimension,dimsf,NULL);

      count[0]=nySub;
      count[1]=nxSub;
      count[2]=nzSub;
      offset[1]=0;
      offset[0]=D->minYSub-D->minYDomain;
      offset[2]=D->minZSub-D->minZDomain;
      memspace=H5Screate_simple(D->dimension,count,NULL);

      field = (float *)malloc(nxSub*nySub*nzSub*sizeof(float ));

      dset_id=H5Dcreate(file_id,"Ex",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      subfilespace=H5Dget_space(dset_id);
      H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
      start=0;
      for(j=jstart; j<jend; j++)
        for(i=istart; i<iend; i++)
        {
          for(k=kstart; k<kend; k++)
            field[start+k-kstart]=D->Ex[i][j][k];
          start+=nzSub;
        }
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
      H5Pclose(plist_id);
      H5Sclose(subfilespace);
      H5Dclose(dset_id);
          
      dset_id=H5Dcreate(file_id,"Ey",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      subfilespace=H5Dget_space(dset_id);
      H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
      start=0;
      for(j=jstart; j<jend; j++)
        for(i=istart; i<iend; i++)
        {
          for(k=kstart; k<kend; k++)
            field[start+k-kstart]=D->Pr[i][j][k]+D->Pl[i][j][k];
          start+=nzSub;
        }
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
      H5Pclose(plist_id);
      H5Sclose(subfilespace);
      H5Dclose(dset_id);

      dset_id=H5Dcreate(file_id,"Ez",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      subfilespace=H5Dget_space(dset_id);
      H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
      start=0;
      for(j=jstart; j<jend; j++)
        for(i=istart; i<iend; i++)
        {
          for(k=kstart; k<kend; k++)
            field[start+k-kstart]=D->Sr[i][j][k]+D->Sl[i][j][k];
          start+=nzSub;
        }
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
      H5Pclose(plist_id);
      H5Sclose(subfilespace);
      H5Dclose(dset_id);

      dset_id=H5Dcreate(file_id,"Bx",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      subfilespace=H5Dget_space(dset_id);
      H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
      start=0;
      for(j=jstart; j<jend; j++)
        for(i=istart; i<iend; i++)
        {
          for(k=kstart; k<kend; k++)
            field[start+k-kstart]=D->Bx[i][j][k];
          start+=nzSub;
        }
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
      H5Pclose(plist_id);
      H5Sclose(subfilespace);
      H5Dclose(dset_id);

      dset_id=H5Dcreate(file_id,"By",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      subfilespace=H5Dget_space(dset_id);
      H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
      start=0;
      for(j=jstart; j<jend; j++)
        for(i=istart; i<iend; i++)
        {
          for(k=kstart; k<kend; k++)
            field[start+k-kstart]=D->Sl[i][j][k]-D->Sr[i][j][k];
          start+=nzSub;
        }
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
      H5Pclose(plist_id);
      H5Sclose(subfilespace);
      H5Dclose(dset_id);

      dset_id=H5Dcreate(file_id,"Bz",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      subfilespace=H5Dget_space(dset_id);
      H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
      start=0;
      for(j=jstart; j<jend; j++)
        for(i=istart; i<iend; i++)
        {
          for(k=kstart; k<kend; k++)
            field[start+k-kstart]=D->Pr[i][j][k]-D->Pl[i][j][k];
          start+=nzSub;
        }
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
      H5Pclose(plist_id);
      H5Sclose(subfilespace);
      H5Dclose(dset_id);

      H5Sclose(memspace);
      H5Sclose(filespace);

      free(field);
    break;
*/
  }		//End of switch(dimension....)
//    H5Fclose(file_id);

//  for(i=0; i<nSpecies; i++)
//    free(cntOffSet[i]);
//  free(cntOffSet);
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

void restoreFloatData(char *fileName,char *dataName,float *data,int totalCnt,int dataNum,int offSet)
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
  status = H5Dread(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,data);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restoreField2D(float ***data,char *fileName,char *dataName,int nx,int ny,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int *offSet)
{
  int i,j,k,start;
  float *field;
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

  field = (float *)malloc(nxSub*nySub*sizeof(float ));

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
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


