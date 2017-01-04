#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void solveDensity2D(Domain *D,int s,float coef)
{
  int i,j,k,ii,jj,istart,iend,jstart,jend;
  float Wx[4],Wy[4];
  float x,x1,x2,x3,x4,y,y1,y2,y3,y4;
  Particle ***particle;
  particle=D->particle;
  ptclList *p;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;

  k=0;
  for(i=0; i<iend+3; i++)
    for(j=0; j<jend+3; j++)
      D->Rho[i][j][k]=0.0;

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
    {
      p=particle[i][j][k].head[s]->pt;
      while(p)
      {
        x=p->x;
        x1=1+x;
        x2=x;
        x3=1-x;
        x4=2-x;
        Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
        Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
        Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
        Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
        y=p->y;
        y1=1+y;
        y2=y;
        y3=1-y;
        y4=2-y;
        Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
        Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
        Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
        Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;
        for(ii=0; ii<4; ii++)
          for(jj=0; jj<4; jj++)
            D->Rho[i-1+ii][j-1+jj][k]+=Wx[ii]*Wy[jj]*coef;
        p=p->next;
      }
    }
}   

void solveDensity3D(Domain *D,int s,float coef)
{
  int i,j,k,ii,jj,kk,istart,iend,jstart,jend,kstart,kend;
  float Wx[4],Wy[4],Wz[4];
  float x,x1,x2,x3,x4,y,y1,y2,y3,y4,z,z1,z2,z3,z4;
  Particle ***particle;
  particle=D->particle;
  ptclList *p;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;

  for(i=0; i<iend+3; i++)
    for(j=0; j<jend+3; j++)
      for(k=0; k<kend+3; k++)
        D->Rho[i][j][k]=0.0;

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      for(k=kstart; k<kend; k++)
      {
        p=particle[i][j][k].head[s]->pt;
        while(p)
        {
          x=p->x;
          x1=1+x;
          x2=x;
          x3=1-x;
          x4=2-x;
          Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
          Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
          Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
          Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
          y=p->y;
          y1=1+y;
          y2=y;
          y3=1-y;
          y4=2-y;
          Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
          Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
          Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
          Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;
          z=p->z;
          z1=1+z;
          z2=z;
          z3=1-z;
          z4=2-z;
          Wz[0]=(2-z1)*(2-z1)*(2-z1)/6.0;
          Wz[1]=(4-6*z2*z2+3*z2*z2*z2)/6.0;
          Wz[2]=(4-6*z3*z3+3*z3*z3*z3)/6.0;
          Wz[3]=(2-z4)*(2-z4)*(2-z4)/6.0;
          for(ii=0; ii<4; ii++)
            for(jj=0; jj<4; jj++)
              for(kk=0; kk<4; kk++)
                D->Rho[i-1+ii][j-1+jj][k-1+kk]+=Wx[ii]*Wy[jj]*Wz[kk]*coef;
          p=p->next;
        }
      }
}  
 
void saveDensityHDF(Domain D,int iteration,int nSpecies)
{
  void density_xdmf();
  void saveDensityCoord_HDF();
  void saveDensity_HDF();

  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  
  saveDensity_HDF(D,iteration);
  MPI_Barrier(MPI_COMM_WORLD);
  saveDensityCoord_HDF(&D,iteration);
  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0)
    density_xdmf(D.dimension,iteration,D.nx,D.ny,D.nz,nSpecies);
  MPI_Barrier(MPI_COMM_WORLD);
}

// The number of cells in the X, Y dimensions
void density_xdmf(int dimension,int iteration,int nx,int ny,int nz,int nSpecies)
{
    FILE *xmf = 0;
    char name[100];
    int s; 
    //
     // Open the file and write the XML description of the mesh..
     //
    sprintf(name,"density%d.xmf",iteration);
    xmf = fopen(name,"w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");

    sprintf(name,"density");    
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
      for(s=0; s<nSpecies; s++)
      {
        fprintf(xmf, "     <Attribute Name=\"%d\" AttributeType=\"Scalar\" Center=\"Node\">\n",s);
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx,ny);
        fprintf(xmf, "        %s%d.h5:/%d\n",name,iteration,s);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
    case 3 :
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
      for(s=0; s<nSpecies; s++)
      {
        fprintf(xmf, "     <Attribute Name=\"%d\" AttributeType=\"Scalar\" Center=\"Node\">\n",s);
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nz,nx,ny);
        fprintf(xmf, "        %s%d.h5:/%d\n",name,iteration,s);
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

void saveDensityCoord_HDF(Domain *D,int iteration)
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
    sprintf(name,"density%d.h5",iteration);
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
  else ;
}

void saveDensity_HDF(Domain D,int iteration)
{
    int i,j,k,s;
    char name[100],dataName[100];
    LoadList *LL;
//    void solveDensity2D();
//    void solveDensity3D();
    int myrank, nTasks;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id;
    herr_t status;
    int offset[3];
    float *rho0;
    float charge;

    sprintf(name,"density%d.h5",iteration);
//    plist_id=H5Pcreate(H5P_FILE_ACCESS);
//    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//    H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else 	;
    s=0;
    rho0 = (float *)malloc((D.nSpecies)*sizeof(float ));
    LL=D.loadList;
    while(LL->next)
    {
      if(LL->charge<0)	charge=-1.0*LL->charge;
      else 	        charge=LL->charge;
      rho0[s]=charge*LL->density/LL->numberInCell;
      LL=LL->next;
      s++;
    }

    switch(D.dimension) {
    //2D
    case 2:
      offset[0]=0;
      offset[1]=D.minYSub-D.minYDomain;
      for(s=0; s<D.nSpecies; s++)
      {
        //solve density
        solveDensity2D(&D,s,rho0[s]);
        if(D.M>1)
        {
          MPI_Density_Yminus(&D,D.nx+5,1,3);
          MPI_Density_Yplus(&D,D.nx+5,1,3);
        }
        else 	;
        sprintf(dataName,"%d",s);
        saveFieldComp_2D(D.Rho,name,dataName,D.nx,D.ny,D.nxSub,D.nySub,D.istart,D.iend,D.jstart,D.jend,offset);        
        MPI_Barrier(MPI_COMM_WORLD);
        if(myrank==0)
        {
          saveIntMeta(name,"nx",&(D.nx));
          saveIntMeta(name,"ny",&(D.ny));
        }
        else	;
        MPI_Barrier(MPI_COMM_WORLD);
      }
      break;
/*
    //3D
    case 3 :
      dimsf[0]=D.ny;      
      dimsf[1]=D.nx;     
      dimsf[2]=D.nz;     
      nx=D.nx;
      ny=D.ny;
      nz=D.nz;
      filespace=H5Screate_simple(D.dimension,dimsf,NULL);

      count[0]=nySub;
      count[1]=nxSub;
      count[2]=nzSub;
      offset[1]=0;
      offset[0]=D.minYSub-D.minYDomain;
      offset[2]=D.minZSub-D.minZDomain;
      memspace=H5Screate_simple(D.dimension,count,NULL);

      totalNum=nxSub*nySub*nzSub;
      density = (double *)malloc(totalNum*sizeof(double ));
//here
      for(s=0; s<D.nSpecies; s++)
      {
        sprintf(name,"%d",s);
        dset_id=H5Dcreate(file_id,name,H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        subfilespace=H5Dget_space(dset_id);
        H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
        //solve density
        solveDensity3D(&D,s,rho0[s]);

        if(D.M>1)
        {
          MPI_Density_Yminus(&D,nx+5,nzSub+5,3);
          MPI_Density_Yplus(&D,nx+5,nzSub+5,3);
        }
        if(D.N>1)
        {
          MPI_Density_Zminus(&D,nx+5,nySub+5,3);
          MPI_Density_Zplus(&D,nx+5,nySub+5,3);
        }

        //asign density data
        start=0;
        for(j=jstart; j<jend; j++)
          for(i=istart; i<iend; i++)
          {
            for(k=kstart; k<kend; k++)
              density[start+k-kstart]=D.Rho[i][j][k];
            start+=nzSub;
          }
        plist_id=H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
        status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,density);
        H5Pclose(plist_id);
        H5Sclose(subfilespace);
        H5Dclose(dset_id);
        MPI_Barrier(MPI_COMM_WORLD);
      }  

      H5Sclose(memspace);
      H5Sclose(filespace);
      free(density);
      break;
*/
    }		//End of switch(dimension....)
    free(rho0);
//    status=H5Fclose(file_id);
}

