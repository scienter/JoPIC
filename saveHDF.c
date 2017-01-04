#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include "hdf5.h"
 
// The number of cells in the X, Y dimensions
 
void field_xdmf(int dimension,int iteration,int nx,int ny,int nz)
{
    FILE *xmf = 0;
    char filename[100];
 
    //
     // Open the file and write the XML description of the mesh..
     //
    sprintf(filename,"field%d.xmf",iteration);
    xmf = fopen(filename,"w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
    
    switch (dimension)  {
    case 2 :
      fprintf(xmf, "     <Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n", ny,nx);
      fprintf(xmf, "     <Geometry GeometryType=\"VXVY\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        field%d.h5:/X\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",ny);
      fprintf(xmf, "        field%d.h5:/Y\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");
      fprintf(xmf, "     <Attribute Name=\"Ex\" AttributeType=\"Scalar\" Center=\"Node\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",ny,nx);
      fprintf(xmf, "        field%d.h5:/Ex\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "     <Attribute Name=\"Ey\" AttributeType=\"Scalar\" Center=\"Node\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",ny,nx);
      fprintf(xmf, "        field%d.h5:/Ey\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "     <Attribute Name=\"Ez\" AttributeType=\"Scalar\" Center=\"Node\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",ny,nx);
      fprintf(xmf, "        field%d.h5:/Ez\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "     <Attribute Name=\"Bx\" AttributeType=\"Scalar\" Center=\"Node\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",ny,nx);
      fprintf(xmf, "        field%d.h5:/Bx\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "     <Attribute Name=\"By\" AttributeType=\"Scalar\" Center=\"Node\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",ny,nx);
      fprintf(xmf, "        field%d.h5:/By\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "     <Attribute Name=\"Bz\" AttributeType=\"Scalar\" Center=\"Node\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",ny,nx);
      fprintf(xmf, "        field%d.h5:/Bz\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
    case 3 :
      fprintf(xmf, "     <Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"%d %d %d\"/>\n",nz,ny,nx);
      fprintf(xmf, "     <Geometry GeometryType=\"VXVYVZ\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        field%d.h5:/X\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",ny);
      fprintf(xmf, "        field%d.h5:/Y\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",nz);
      fprintf(xmf, "        field%d.h5:/Z\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");
      fprintf(xmf, "     <Attribute Name=\"Ex\" AttributeType=\"Scalar\" Center=\"Node\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",nx,ny,nz);
      fprintf(xmf, "        field%d.h5:/Ex\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "     <Attribute Name=\"Ey\" AttributeType=\"Scalar\" Center=\"Node\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",nx,ny,nz);
      fprintf(xmf, "        field%d.h5:/Ey\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "     <Attribute Name=\"Ez\" AttributeType=\"Scalar\" Center=\"Node\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",nx,ny,nz);
      fprintf(xmf, "        field%d.h5:/Ez\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "     <Attribute Name=\"Bx\" AttributeType=\"Scalar\" Center=\"Node\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",nx,ny,nz);
      fprintf(xmf, "        field%d.h5:/Bx\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "     <Attribute Name=\"By\" AttributeType=\"Scalar\" Center=\"Node\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",nx,ny,nz);
      fprintf(xmf, "        field%d.h5:/By\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "     <Attribute Name=\"Bz\" AttributeType=\"Scalar\" Center=\"Node\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",nx,ny,nz);
      fprintf(xmf, "        field%d.h5:/Bz\n",iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
    }
}

void saveField_HDF(Domain *D,int iteration)
{
    int ii,i,j,k,start,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    char name[100];
    double x,y,z,Ex,Ey,Ez,Bx,By,Bz,factor;
    double *field,*xtic,*ytic,*ztic;
    FILE *out,*out1;
    int myrank, nTasks;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[3],count[3],offset[3],dimy[1],dimx[1],dimz[1];
    const char *coorName[] = {"/Y","/X","/Z"};  

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
      dimsf[0]=D->ny;      
      dimsf[1]=D->nx;     
      nx=D->nx;
      ny=D->ny;

      filespace=H5Screate_simple(D->dimension,dimsf,NULL);

      count[0]=nySub;
      count[1]=nxSub;
      offset[1]=0;
      offset[0]=D->minYSub-D->minYDomain;
      memspace=H5Screate_simple(D->dimension,count,NULL);
      field = (double *)malloc(nxSub*nySub*sizeof(double ));

      for(ii=0; ii<6; ii++)
      {
        switch(ii)  {
        case 0:
          dset_id=H5Dcreate(file_id,"Ex",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
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
          status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
          H5Pclose(plist_id);
          H5Sclose(subfilespace);
          H5Dclose(dset_id);
          break;
        case 1:
          dset_id=H5Dcreate(file_id,"Ey",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
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
          status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
          H5Pclose(plist_id);
          H5Sclose(subfilespace);
          H5Dclose(dset_id);
          break;
        case 2:
          dset_id=H5Dcreate(file_id,"Ez",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
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
          status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
          H5Pclose(plist_id);
          H5Sclose(subfilespace);
          H5Dclose(dset_id);
          break;
        case 3:
          dset_id=H5Dcreate(file_id,"Bx",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
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
          status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
          H5Pclose(plist_id);
          H5Sclose(subfilespace);
          H5Dclose(dset_id);
          break;
        case 4:
          dset_id=H5Dcreate(file_id,"By",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
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
          status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
          H5Pclose(plist_id);
          H5Sclose(subfilespace);
          H5Dclose(dset_id);
          break;
        case 5:
          dset_id=H5Dcreate(file_id,"Bz",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
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
          status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
          H5Pclose(plist_id);
          H5Sclose(subfilespace);
          H5Dclose(dset_id);
          break;
        }
      }
      H5Sclose(filespace);

      dimx[0]=nx;
      dimy[0]=ny;
      xtic=(double *)malloc(nx*sizeof(double));
      for(i=0;i<nx;i++)
        xtic[i]=(i+D->minXSub)*D->lambda*D->dx;
      ytic=(double *)malloc(ny*sizeof(double));
      for(i=0;i<ny;i++)
        ytic[i]=(i+D->minYDomain)*D->lambda*D->dy;
      for(ii=0; ii<2; ii++)
      {
        if(ii==0)
          filespace=H5Screate_simple(1,dimy,NULL);
        else if(ii==1)
          filespace=H5Screate_simple(1,dimx,NULL);
        dset_id=H5Dcreate(file_id,coorName[ii],H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,ii==0 ? ytic : xtic);
        H5Dclose(dset_id);
        H5Sclose(filespace);
      }

      free(xtic);
      free(ytic);
      free(field);
      H5Sclose(memspace);
    break;

    //3D
    case (Split-1)*3+3:
      dimsf[0]=D->nz;      
      dimsf[1]=D->ny;     
      dimsf[2]=D->nx;     
      nx=D->nx;
      ny=D->ny;
      nz=D->nz;

      filespace=H5Screate_simple(D->dimension,dimsf,NULL);

      count[0]=nzSub;
      count[1]=nySub;
      count[2]=nxSub;
      offset[2]=0;
      offset[1]=D->minYSub-D->minYDomain;
      offset[0]=D->minZSub-D->minZDomain;
      memspace=H5Screate_simple(D->dimension,count,NULL);
      field = (double *)malloc(nxSub*nySub*nzSub*sizeof(double ));

      dset_id=H5Dcreate(file_id,"Ex",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      subfilespace=H5Dget_space(dset_id);
      H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
      start=0;
      for(k=kstart; k<kend; k++)
        for(j=jstart; j<jend; j++)
        {
          for(i=istart; i<iend; i++)
            field[start+i-istart]=D->Ex[i][j][k];
          start+=nxSub;
        }
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
      H5Pclose(plist_id);
      H5Sclose(subfilespace);
      H5Dclose(dset_id);
          
      dset_id=H5Dcreate(file_id,"Ey",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      subfilespace=H5Dget_space(dset_id);
      H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
      start=0;
      for(k=kstart; k<kend; k++)
        for(j=jstart; j<jend; j++)
        {
          for(i=istart; i<iend; i++)
            field[start+i-istart]=D->Pr[i][j][k]+D->Pl[i][j][k];
          start+=nxSub;
        }
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
      H5Pclose(plist_id);
      H5Sclose(subfilespace);
      H5Dclose(dset_id);

      dset_id=H5Dcreate(file_id,"Ez",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      subfilespace=H5Dget_space(dset_id);
      H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
      start=0;
      for(k=kstart; k<kend; k++)
        for(j=jstart; j<jend; j++)
        {
          for(i=istart; i<iend; i++)
            field[start+i-istart]=D->Sr[i][j][k]+D->Sl[i][j][k];
          start+=nxSub;
        }
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
      H5Pclose(plist_id);
      H5Sclose(subfilespace);
      H5Dclose(dset_id);

      dset_id=H5Dcreate(file_id,"Bx",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      subfilespace=H5Dget_space(dset_id);
      H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
      start=0;
      for(k=kstart; k<kend; k++)
        for(j=jstart; j<jend; j++)
        {
          for(i=istart; i<iend; i++)
            field[start+i-istart]=D->Bx[i][j][k];
          start+=nxSub;
        }
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
      H5Pclose(plist_id);
      H5Sclose(subfilespace);
      H5Dclose(dset_id);

      dset_id=H5Dcreate(file_id,"By",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      subfilespace=H5Dget_space(dset_id);
      H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
      start=0;
      for(k=kstart; k<kend; k++)
        for(j=jstart; j<jend; j++)
        {
          for(i=istart; i<iend; i++)
            field[start+i-istart]=D->Sl[i][j][k]-D->Sr[i][j][k];
          start+=nxSub;
        }
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
      H5Pclose(plist_id);
      H5Sclose(subfilespace);
      H5Dclose(dset_id);

      dset_id=H5Dcreate(file_id,"Bz",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      subfilespace=H5Dget_space(dset_id);
      H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
      start=0;
      for(k=kstart; k<kend; k++)
        for(j=jstart; j<jend; j++)
        {
          for(i=istart; i<iend; i++)
            field[start+i-istart]=D->Pr[i][j][k]-D->Pl[i][j][k];
          start+=nxSub;
        }
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
      H5Pclose(plist_id);
      H5Sclose(subfilespace);
      H5Dclose(dset_id);
      H5Sclose(filespace);


      dimx[0]=nx;
      dimy[0]=ny;
      dimz[0]=nz;
      xtic=(double *)malloc(nx*sizeof(double));
      for(i=0;i<nx;i++)
        xtic[i]=(i+D->minXSub)*D->lambda*D->dx;
      ytic=(double *)malloc(ny*sizeof(double));
      for(i=0;i<ny;i++)
        ytic[i]=(i+D->minYDomain)*D->lambda*D->dy;
      ztic=(double *)malloc(nz*sizeof(double));
      for(i=0;i<nz;i++)
        ztic[i]=(i+D->minZDomain)*D->lambda*D->dz;
      filespace=H5Screate_simple(1,dimy,NULL);
      dset_id=H5Dcreate(file_id,"/Y",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,ytic);
      H5Dclose(dset_id);
      H5Sclose(filespace);
      filespace=H5Screate_simple(1,dimx,NULL);
      dset_id=H5Dcreate(file_id,"/X",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,xtic);
      H5Dclose(dset_id);
      H5Sclose(filespace);
      filespace=H5Screate_simple(1,dimz,NULL);
      dset_id=H5Dcreate(file_id,"/Z",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,ztic);
      H5Dclose(dset_id);
      H5Sclose(filespace);
      free(xtic);
      free(ytic);
      free(ztic);

      free(field);
      H5Sclose(memspace);
    break;
    }		//End of switch(dimension....)

    status=H5Fclose(file_id);
}

