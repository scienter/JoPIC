#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"
 
void saveRamanHDF(Domain D,int iteration)
{
  void raman_xdmf();
  void saveRamanCoord_HDF();
  void saveRaman_HDF();
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  saveRaman_HDF(&D,iteration);
  MPI_Barrier(MPI_COMM_WORLD);
  saveRamanCoord_HDF(&D,iteration);
  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0)
    raman_xdmf(D.dimension,iteration,D.nx,D.ny,D.nz);
  else ;
  MPI_Barrier(MPI_COMM_WORLD);
}

void raman_xdmf(int dimension,int iteration,int nx,int ny,int nz)
{
    FILE *xmf = 0;
    char name[100];
    const char *Names[] = {"/Pr","/Pl","/Sr","/Sl"};
    int i;
    //
     // Open the file and write the XML description of the mesh..
     //
    sprintf(name,"raman%d.xmf",iteration);
    xmf = fopen(name,"w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
    
    sprintf(name,"raman");
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
      for(i=0; i<4; i++)
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
      for(i=0; i<4; i++)
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

void saveRamanCoord_HDF(Domain *D,int iteration)
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
    sprintf(name,"raman%d.h5",iteration);
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


void saveRaman_HDF(Domain *D,int iteration)
{
    char name[100];
    int myrank, nTasks;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id;
    int offset[3];

    sprintf(name,"raman%d.h5",iteration);
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
      offset[0]=0;
      offset[1]=D->minYSub-D->minYDomain;
      saveFieldComp_2D(D->Pr,name,"/Pr",D->nx,D->ny,D->nxSub,D->nySub,D->istart,D->iend,D->jstart,D->jend,offset);
      saveFieldComp_2D(D->Pl,name,"/Pl",D->nx,D->ny,D->nxSub,D->nySub,D->istart,D->iend,D->jstart,D->jend,offset);
      saveFieldComp_2D(D->Sr,name,"/Sr",D->nx,D->ny,D->nxSub,D->nySub,D->istart,D->iend,D->jstart,D->jend,offset);
      saveFieldComp_2D(D->Sl,name,"/Sl",D->nx,D->ny,D->nxSub,D->nySub,D->istart,D->iend,D->jstart,D->jend,offset);
      break;

    //3D
    case (Split-1)*3+3:
      offset[1]=0;
      offset[0]=D->minYSub-D->minYDomain;
      offset[2]=D->minZSub-D->minZDomain;
      saveFieldComp_3D(D->Pr,name,"/Pr",D->nx,D->ny,D->nz,D->nxSub,D->nySub,D->nzSub,D->istart,D->iend,D->jstart,D->jend,D->kstart,D->kend,offset);
      saveFieldComp_3D(D->Pl,name,"/Pl",D->nx,D->ny,D->nz,D->nxSub,D->nySub,D->nzSub,D->istart,D->iend,D->jstart,D->jend,D->kstart,D->kend,offset);
      saveFieldComp_3D(D->Sr,name,"/Sr",D->nx,D->ny,D->nz,D->nxSub,D->nySub,D->nzSub,D->istart,D->iend,D->jstart,D->jend,D->kstart,D->kend,offset);
      saveFieldComp_3D(D->Sl,name,"/Sl",D->nx,D->ny,D->nz,D->nxSub,D->nySub,D->nzSub,D->istart,D->iend,D->jstart,D->jend,D->kstart,D->kend,offset);
      break;
    }		//End of switch(dimension....)
}

