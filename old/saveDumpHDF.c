#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"
 
void saveDumpHDF(Domain D,int iteration)
{
  void saveMeta_HDF();
  void saveDump_HDF();
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  saveDump_HDF(&D,iteration);
  MPI_Barrier(MPI_COMM_WORLD);
//  saveFieldCoord_HDF(&D,iteration);
//  MPI_Barrier(MPI_COMM_WORLD);
}
/*
void saveFieldCoord_HDF(Domain *D,int iteration)
{
  int ii,i,nx,ny,nz;
  char name[100];
  float *xtic,*ytic,*ztic;
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,tic_id;
  herr_t status;
  hid_t filespace;
  hsize_t dimy[1],dimx[1],dimz[1];
  const char *coorName[] = {"/Y","/X","/Z"};  

  nx=D->nx;
  ny=D->ny;

  if(myrank==0)
  {
    sprintf(name,"field%d.h5",iteration);
    file_id=H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);

    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Split-1)*3+2:
      dimx[0]=nx;
      dimy[0]=ny;
      xtic=(float *)malloc(nx*sizeof(float));
      for(i=0;i<nx;i++)
        xtic[i]=(i+D->minXSub)*D->lambda*D->dx;
      ytic=(float *)malloc(ny*sizeof(float));
      for(i=0;i<ny;i++)
        ytic[i]=(i+D->minYDomain)*D->lambda*D->dy;
      for(ii=0; ii<2; ii++)
      {
        if(ii==0)
          filespace=H5Screate_simple(1,dimy,NULL);
        else if(ii==1)
          filespace=H5Screate_simple(1,dimx,NULL);
        dset_id=H5Dcreate(file_id,coorName[ii],H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,ii==0 ? ytic : xtic);
        H5Dclose(dset_id);
        H5Sclose(filespace);
      }
      free(xtic);
      free(ytic);
      break;
    case (Split-1)*3+3:
      dimx[0]=nx;
      dimy[0]=ny;
      dimz[0]=nz;
      xtic=(float *)malloc(nx*sizeof(float));
      for(i=0;i<nx;i++)
        xtic[i]=(i+D->minXSub)*D->lambda*D->dx;
      ytic=(float *)malloc(ny*sizeof(float));
      for(i=0;i<ny;i++)
        ytic[i]=(i+D->minYDomain)*D->lambda*D->dy;
      ztic=(float *)malloc(nz*sizeof(float));
      for(i=0;i<nz;i++)
        ztic[i]=(i+D->minZDomain)*D->lambda*D->dz;
      filespace=H5Screate_simple(1,dimy,NULL);
      dset_id=H5Dcreate(file_id,"/Y",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,ytic);
      H5Dclose(dset_id);
      H5Sclose(filespace);
      filespace=H5Screate_simple(1,dimx,NULL);
      dset_id=H5Dcreate(file_id,"/X",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,xtic);
      H5Dclose(dset_id);
      H5Sclose(filespace);
      filespace=H5Screate_simple(1,dimz,NULL);
      dset_id=H5Dcreate(file_id,"/Z",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,ztic);
      H5Dclose(dset_id);
      H5Sclose(filespace);

      free(xtic);
      free(ytic);
      free(ztic);
      break;
    }
    H5Fclose(file_id);
  }
  else  ;
}
*/

void saveDump_HDF(Domain *D,int iteration)
{
    int i,istart,iend,jstart,jend,kstart,kend,nx,ny,nz,bias,offSetY;
    int nxSub,nySub,nzSub,*offSetJ;
    char name[100];
    void saveFileField();

    int myrank, nTasks;    
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
   
    metaDim[0]=nTasks;
  
    sprintf(name,"dump%d.h5",iteration);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Split-1)*3+2:
      offSetJ = (int *)malloc(nTasks*sizeof(int ));

      ny=D->ny+5;      
      nx=D->nx+5;     
      if(myrank==0)  {      
        jstart=0;
        nySub+=2;    
        bias=0;
      }
      else if(myrank==nTasks-1)  {
        jend+=3;
        nySub+=3; 
        bias=2;
      }
      else   {
        bias=2; 
      } 
      istart=0;
      iend+=3; 
      nxSub+=5;
      offSetY=D->minYSub-D->minYDomain+bias;
      
      MPI_Gather(&offSetY,1,MPI_INT,offSetJ,1,MPI_INT,0,MPI_COMM_WORLD);
      if(myrank==0)
      {  
        sprintf(name,"dump%d.h5",iteration);
        file_id=H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);
        filespace=H5Screate_simple(1,metaDim,NULL);
        dset_id=H5Dcreate(file_id,"/offSetJ",H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,offSetJ);
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Fclose(file_id);
      }
      else  	;
     

      saveFileField(D->Ex,"/Ex",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->Pr,"/Pr",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->Pl,"/Pl",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->Bx,"/Bx",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->Sr,"/Sr",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->Sl,"/Sl",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->ExC,"/ExC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->PrC,"/PrC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->PlC,"/PlC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->BxC,"/BxC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->SrC,"/SrC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->SlC,"/SlC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->Jx,"/Jx",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->Jy,"/Jy",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->Jz,"/Jz",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->JxOld,"/JxOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->JyOld,"/JyOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);
      saveFileField(D->JzOld,"/JzOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,kstart,kend,iteration,D->dimension,offSetY);

      free(offSetJ);
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

}


void saveFileField(float ***data,char *dataName,int nx,int ny,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int kstart,int kend,int iteration,int dimension,int offSetY)
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

  sprintf(name,"dump%d.h5",iteration);
  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
  H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
  MPI_Barrier(MPI_COMM_WORLD);

  file_id=H5Fopen(name,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[1]=nx;      
  dimsf[0]=ny;     
  filespace=H5Screate_simple(dimension,dimsf,NULL);

  count[0]=nySub;
  count[1]=nxSub;
  offset[1]=0;
  offset[0]=offSetY;
  memspace=H5Screate_simple(dimension,count,NULL);

  field = (float *)malloc(nxSub*nySub*sizeof(float ));

  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  start=0;
  for(j=jstart; j<jend; j++)
  {
    for(i=istart; i<iend; i++)
      field[start+i-istart]=data[i][j][0];
    start+=nxSub;
  }
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);
  MPI_Barrier(MPI_COMM_WORLD);

  H5Sclose(memspace);
  H5Sclose(filespace);
  free(field);
  H5Fclose(file_id);
}
