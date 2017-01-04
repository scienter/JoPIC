#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"
 
void saveFieldHDF(Domain D,int iteration)
{
  void createField();
  void field_xdmf();
  void saveFieldCoord_HDF();
  void saveField_HDF();
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  saveField_HDF(&D,iteration);
  MPI_Barrier(MPI_COMM_WORLD);
  saveFieldCoord_HDF(&D,iteration);
  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0)
    field_xdmf(D.dimension,iteration,D.nx,D.ny,D.nz);
  else ;
  MPI_Barrier(MPI_COMM_WORLD);
}

void field_xdmf(int dimension,int iteration,int nx,int ny,int nz)
{
    FILE *xmf = 0;
    char name[100];
    const char *Names[] = {"/Ex","/Ey","/Ez","/Bx","/By","/Bz"};
    int i;
    //
     // Open the file and write the XML description of the mesh..
     //
    sprintf(name,"field%d.xmf",iteration);
    xmf = fopen(name,"w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
    
    sprintf(name,"field");
    switch (dimension)  {
    case 2 :
      fprintf(xmf, "     <Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n", nx,ny);
      fprintf(xmf, "     <Geometry GeometryType=\"VXVY\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ny);
      fprintf(xmf, "        %s%d.h5:/Y\n",name,iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s%d.h5:/X\n",name,iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");
      for(i=0; i<6; i++)
      {
        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Names[i]);
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx,ny);
        fprintf(xmf, "        %s%d.h5:/%s\n",name,iteration,Names[i]);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
    case 3 :
      //3D
      fprintf(xmf, "     <Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"%d %d %d\"/>\n",ny,nx,nz);
      fprintf(xmf, "     <Geometry GeometryType=\"VXVYVZ\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nz);
      fprintf(xmf, "        %s%d.h5:/Z\n",name,iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s%d.h5:/X\n",name,iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ny);
      fprintf(xmf, "        %s%d.h5:/Y\n",name,iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");   
      for(i=0; i<6; i++)
      {
        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Names[i]);
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nz,nx,ny);
        fprintf(xmf, "        %s%d.h5:/%s\n",name,iteration,Names[i]);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
    }
}

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
  nz=D->nz;

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


void saveField_HDF(Domain *D,int iteration)
{
    int ii,i,j,k,start,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    char name[100],dataName[100];
    float x,y,z,Ex,Ey,Ez,Bx,By,Bz,factor;
    float *field,*xtic,*ytic,*ztic;
    FILE *out,*out1;
    void saveIntMeta();
    int myrank, nTasks;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[3],count[3],dimy[1],dimx[1],dimz[1];
//    hsize_t dimsf[2],count[2],offset[2],dimy[1],dimx[1],dimz[1];

    int did,size,nxSub,nySub,nzSub;
    int *offset;
    void saveFieldComp_2D();
    void saveFieldComp2_2D();
    void saveFieldComp_3D();
    void saveFieldComp2_3D();

    offset=(int *)malloc(3*sizeof(int));
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;
    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;

    sprintf(name,"field%d.h5",iteration);
//    plist_id=H5Pcreate(H5P_FILE_ACCESS);
//    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//    H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else	;
    MPI_Barrier(MPI_COMM_WORLD);

    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Split-1)*3+2:
      if(myrank==0)
      {
        sprintf(dataName,"/nx");
        saveIntMeta(name,dataName,&(D->nx));
        sprintf(dataName,"/ny");
        saveIntMeta(name,dataName,&(D->ny));
      }
      else    ;

      offset[0]=0;
      offset[1]=D->minYSub-D->minYDomain;
      saveFieldComp_2D(D->Ex,name,"/Ex",D->nx,D->ny,D->nxSub,D->nySub,D->istart,D->iend,D->jstart,D->jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D->Bx,name,"/Bx",D->nx,D->ny,D->nxSub,D->nySub,D->istart,D->iend,D->jstart,D->jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp2_2D(D->Pr,D->Pl,name,"/Ey",D->nx,D->ny,D->nxSub,D->nySub,D->istart,D->iend,D->jstart,D->jend,offset,0);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp2_2D(D->Sr,D->Sl,name,"/Ez",D->nx,D->ny,D->nxSub,D->nySub,D->istart,D->iend,D->jstart,D->jend,offset,0);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp2_2D(D->Sl,D->Sr,name,"/By",D->nx,D->ny,D->nxSub,D->nySub,D->istart,D->iend,D->jstart,D->jend,offset,1);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp2_2D(D->Pr,D->Pl,name,"/Bz",D->nx,D->ny,D->nxSub,D->nySub,D->istart,D->iend,D->jstart,D->jend,offset,1);
      MPI_Barrier(MPI_COMM_WORLD);
      break;

    //3D
    case (Split-1)*3+3:
      if(myrank==0)
      {
        sprintf(dataName,"/nx");
        saveIntMeta(name,dataName,&(D->nx));
        sprintf(dataName,"/ny");
        saveIntMeta(name,dataName,&(D->ny));
        sprintf(dataName,"/nz");
        saveIntMeta(name,dataName,&(D->nz));
      }
      else    ;

      offset[1]=0;
      offset[0]=D->minYSub-D->minYDomain;
      offset[2]=D->minZSub-D->minZDomain;
      saveFieldComp_3D(D->Ex,name,"/Ex",D->nx,D->ny,D->nz,D->nxSub,D->nySub,D->nzSub,D->istart,D->iend,D->jstart,D->jend,D->kstart,D->kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D->Bx,name,"/Bx",D->nx,D->ny,D->nz,D->nxSub,D->nySub,D->nzSub,D->istart,D->iend,D->jstart,D->jend,D->kstart,D->kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp2_3D(D->Pr,D->Pl,name,"/Ey",D->nx,D->ny,D->nz,D->nxSub,D->nySub,D->nzSub,D->istart,D->iend,D->jstart,D->jend,D->kstart,D->kend,offset,0);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp2_3D(D->Sr,D->Sl,name,"/Ez",D->nx,D->ny,D->nz,D->nxSub,D->nySub,D->nzSub,D->istart,D->iend,D->jstart,D->jend,D->kstart,D->kend,offset,0);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp2_3D(D->Sl,D->Sr,name,"/By",D->nx,D->ny,D->nz,D->nxSub,D->nySub,D->nzSub,D->istart,D->iend,D->jstart,D->jend,D->kstart,D->kend,offset,1);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp2_3D(D->Pr,D->Pl,name,"/Bz",D->nx,D->ny,D->nz,D->nxSub,D->nySub,D->nzSub,D->istart,D->iend,D->jstart,D->jend,D->kstart,D->kend,offset,1);
      MPI_Barrier(MPI_COMM_WORLD);
    break;

    }		//End of switch(dimension....)
  free(offset);
}


void saveFieldComp_2D(double ***data,char *fileName,char *dataName,int nx,int ny,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int *offSet)
{
    int ii,i,j,k,start;
    float *field;
    int myrank, nTasks;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[2],count[2],offset[2];

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//    H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//    MPI_Barrier(MPI_COMM_WORLD);

    file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);

    dimsf[1]=ny;      
    dimsf[0]=nx;     
    filespace=H5Screate_simple(2,dimsf,NULL);

    count[1]=nySub;
    count[0]=nxSub;
    offset[1]=offSet[1];
    offset[0]=offSet[0];
    memspace=H5Screate_simple(2,count,NULL);

    field = (float *)malloc(nxSub*nySub*sizeof(float ));

    dset_id=H5Dcreate(file_id,dataName,H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    subfilespace=H5Dget_space(dset_id);
    H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
    start=0;
    for(i=istart; i<iend; i++)
    {
      for(j=jstart; j<jend; j++)
        field[start+j-jstart]=data[i][j][0];
      start+=nySub;
    }
    plist_id=H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
    H5Pclose(plist_id);
    H5Sclose(subfilespace);
    H5Dclose(dset_id);
      
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Fclose(file_id);
    free(field);
}

void saveFieldComp2_2D(double ***data1,double ***data2,char *fileName,char *dataName,int nx,int ny,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int *offSet,int mode)
{
    int ii,i,j,k,start;
    float x,y,z,Ex,Ey,Ez,Bx,By,Bz,factor;
    float *field,*xtic,*ytic,*ztic;
    FILE *out,*out1;
    int myrank, nTasks;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[2],count[2],offset[2];

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//    H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//    MPI_Barrier(MPI_COMM_WORLD);

    file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);

    dimsf[1]=ny;      
    dimsf[0]=nx;     
    filespace=H5Screate_simple(2,dimsf,NULL);

    count[1]=nySub;
    count[0]=nxSub;
    offset[1]=offSet[1];
    offset[0]=offSet[0];
    memspace=H5Screate_simple(2,count,NULL);

    field = (float *)malloc(nxSub*nySub*sizeof(float ));

    dset_id=H5Dcreate(file_id,dataName,H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    subfilespace=H5Dget_space(dset_id);
    H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
    start=0;
    if(mode==0)
    {
      for(i=istart; i<iend; i++)
      {
        for(j=jstart; j<jend; j++)
          field[start+j-jstart]=data1[i][j][0]+data2[i][j][0];
        start+=nySub;
      }
    }
    else if(mode==1)
    {
      for(i=istart; i<iend; i++)
      {
        for(j=jstart; j<jend; j++)
          field[start+j-jstart]=data1[i][j][0]-data2[i][j][0];
        start+=nySub;
      }
    }
    else 	;
    plist_id=H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
    H5Pclose(plist_id);
    H5Sclose(subfilespace);
    H5Dclose(dset_id);
      
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Fclose(file_id);
    free(field);
}

void saveFieldComp_3D(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet)
{
    int ii,i,j,k,start;
    float x,y,z,Ex,Ey,Ez,Bx,By,Bz,factor;
    float *field,*xtic,*ytic,*ztic;
    FILE *out,*out1;
    int myrank, nTasks;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[3],count[3],offset[3];

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//    H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//    MPI_Barrier(MPI_COMM_WORLD);

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

    field = (float *)malloc(nxSub*nySub*nzSub*sizeof(float ));

    dset_id=H5Dcreate(file_id,dataName,H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    subfilespace=H5Dget_space(dset_id);
    H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
    start=0;
    for(j=jstart; j<jend; j++)
      for(i=istart; i<iend; i++)
      {
        for(k=kstart; k<kend; k++)
          field[start+k-kstart]=data[i][j][k];
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
    H5Fclose(file_id);
    free(field);
}

void saveFieldComp2_3D(double ***data1,double ***data2,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet,int mode)
{
    int ii,i,j,k,start;
    float x,y,z,Ex,Ey,Ez,Bx,By,Bz,factor;
    float *field,*xtic,*ytic,*ztic;
    FILE *out,*out1;
    int myrank, nTasks;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[3],count[3],offset[3];

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//    H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//    MPI_Barrier(MPI_COMM_WORLD);

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

    field = (float *)malloc(nxSub*nySub*nzSub*sizeof(float ));

    dset_id=H5Dcreate(file_id,dataName,H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    subfilespace=H5Dget_space(dset_id);
    H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
    start=0;
    if(mode==0)
    {
      for(j=jstart; j<jend; j++)
        for(i=istart; i<iend; i++)
        {
          for(k=kstart; k<kend; k++)
            field[start+k-kstart]=data1[i][j][k]+data2[i][j][k];
          start+=nzSub;
        }
    }
    else if(mode==1)
    {
      for(j=jstart; j<jend; j++)
        for(i=istart; i<iend; i++)
        {
          for(k=kstart; k<kend; k++)
            field[start+k-kstart]=data1[i][j][k]-data2[i][j][k];
          start+=nzSub;
        }
    }
    else 	;
    plist_id=H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,field);
    H5Pclose(plist_id);
    H5Sclose(subfilespace);
    H5Dclose(dset_id);
      
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Fclose(file_id);
    free(field);
}

void saveIntMeta(char *fileName,char *dataName,int *data)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=1;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dcreate(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}
