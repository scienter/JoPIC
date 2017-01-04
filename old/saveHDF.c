#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void saveField_HDF(Domain *D,int iteration)
{
    int ii,i,j,k,start,istart,iend,jstart,jend,kstart,kend,nx,ny;
    char name[100];
    float x,y,z,Ex,Ey,Ez,Bx,By,Bz,factor;
    float *field,*xtic,*ytic;
    FILE *out,*out1;
    int myrank, nTasks;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[3],count[3],offset[3];

    const char *Names[] = {"/Ex","/Ey","/Ez","/Bx","/By","/Bz"};
    int did,size,nxSub,nySub,nzSub;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;
    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
    sprintf(name,"field%d.h5",iteration);
    file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
    H5Pclose(plist_id);

    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Split-1)*3+2:
      if(D->M>1)
      {
        dimsf[0]=D->ny;      
        dimsf[1]=D->nx;     
      }
      else
      {
        dimsf[0]=nySub;      
        dimsf[1]=nxSub;
      }     
      nx=D->nx;
      ny=D->ny;

      filespace=H5Screate_simple(D->dimension,dimsf,NULL);

      count[0]=nySub;
      count[1]=nxSub;
      offset[1]=0;
      offset[0]=D->minYSub-D->minYDomain;
      memspace=H5Screate_simple(D->dimension,count,NULL);
      field = (float *)malloc(nxSub*nySub*sizeof(float ));

      for(ii=0; ii<6; ii++)
      {
        switch(ii)  {
        case 0:
          dset_id=H5Dcreate(file_id,"Ex",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          subfilespace=H5Dget_space(dset_id);
          H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
          start=0;
          for(j=jstart; j<jend; j++)
          {
            for(i=istart; i<iend; i++)
              field[start+i-istart]=D->Ex[i][j][0];
            start+=nxSub;
          }
          plist_id=H5Pcreate(H5P_DATASET_XFER);
          H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
          status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
          H5Pclose(plist_id);
          H5Sclose(subfilespace);
          H5Dclose(dset_id);
          break;
        case 1:
          dset_id=H5Dcreate(file_id,"Ey",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          subfilespace=H5Dget_space(dset_id);
          H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
          start=0;
          for(j=jstart; j<jend; j++)
          {
            for(i=istart; i<iend; i++)
              field[start+i-istart]=D->Pr[i][j][0]+D->Pl[i][j][0];
            start+=nxSub;
          }
          plist_id=H5Pcreate(H5P_DATASET_XFER);
          H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
          status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
          H5Pclose(plist_id);
          H5Sclose(subfilespace);
          H5Dclose(dset_id);
          break;
        case 2:
          dset_id=H5Dcreate(file_id,"Ez",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          subfilespace=H5Dget_space(dset_id);
          H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
          start=0;
          for(j=jstart; j<jend; j++)
          {
            for(i=istart; i<iend; i++)
              field[start+i-istart]=D->Sr[i][j][0]+D->Sl[i][j][0];
            start+=nxSub;
          }
          plist_id=H5Pcreate(H5P_DATASET_XFER);
          H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
          status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
          H5Pclose(plist_id);
          H5Sclose(subfilespace);
          H5Dclose(dset_id);
          break;
        case 3:
          dset_id=H5Dcreate(file_id,"Bx",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          subfilespace=H5Dget_space(dset_id);
          H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
          start=0;
          for(j=jstart; j<jend; j++)
          {
            for(i=istart; i<iend; i++)
              field[start+i-istart]=D->Bx[i][j][0];
            start+=nxSub;
          }
          plist_id=H5Pcreate(H5P_DATASET_XFER);
          H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
          status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
          H5Pclose(plist_id);
          H5Sclose(subfilespace);
          H5Dclose(dset_id);
          break;
        case 4:
          dset_id=H5Dcreate(file_id,"By",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          subfilespace=H5Dget_space(dset_id);
          H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
          start=0;
          for(j=jstart; j<jend; j++)
          {
            for(i=istart; i<iend; i++)
              field[start+i-istart]=D->Sl[i][j][0]-D->Sr[i][j][0];
            start+=nxSub;
          }
          plist_id=H5Pcreate(H5P_DATASET_XFER);
          H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
          status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
          H5Pclose(plist_id);
          H5Sclose(subfilespace);
          H5Dclose(dset_id);
          break;
        case 5:
          dset_id=H5Dcreate(file_id,"Bz",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          subfilespace=H5Dget_space(dset_id);
          H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
          start=0;
          for(j=jstart; j<jend; j++)
          {
            for(i=istart; i<iend; i++)
              field[start+i-istart]=D->Pr[i][j][0]-D->Pl[i][j][0];
            start+=nxSub;
          }
          plist_id=H5Pcreate(H5P_DATASET_XFER);
          H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
          status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
          H5Pclose(plist_id);
          H5Sclose(subfilespace);
          break;
        }
      }

      hsize_t x_dim[1]={nx};
      xtic=(float *)malloc(nx*sizeof(float));
      for(i=0;i<nx;i++)
        xtic[i]=(i+D->minXSub)*D->lambda*D->dx;
      H5LTmake_dataset_float(file_id,"x",1,x_dim,xtic);
      tic_id=H5Dopen2(file_id,"x",H5P_DEFAULT);
      H5DSattach_scale(dset_id,tic_id,1);	//1:index of x      
      H5Dclose(tic_id);
      
      hsize_t y_dim[1]={ny};
      ytic=(float *)malloc(ny*sizeof(float));
      for(i=0;i<ny;i++)
        ytic[i]=(i+D->minYDomain)*D->lambda*D->dy;
      H5LTmake_dataset_float(file_id,"y",1,y_dim,ytic);
      tic_id=H5Dopen2(file_id,"y",H5P_DEFAULT);
      H5DSattach_scale(dset_id,tic_id,0);	//0:index of y      
      H5Dclose(tic_id);

      H5Dclose(dset_id);
      H5Sclose(filespace);
      free(xtic);
      free(ytic);
      H5Sclose(memspace);
      free(field);
    break;
    }		//End of switch(dimension....)

    status=H5Fclose(file_id);
}

