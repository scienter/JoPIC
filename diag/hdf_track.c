#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "mpi.h"
#include "math.h"

void restoreData(char *fileName,char *dataName,int totalCnt,int cntSub,int start,double *data,int startC,int columns);
void restore1Data(char *fileName,char *dataName,int totalCnt,int cntSub,int start,double *data,int startC);
void calSub(int totalCnt,int *cntSub,int *start);
void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
double whatMass(int atomNum);

void main(int argc, char *argv[])
{
   int mode,species,step,atomNum,core;
   int i,n,totalCnt,cntSub,start,column,sharePNum,index,dataNum=9;
   int *recv;
   double minPz,maxPz,minZ,maxZ,aveGam,cnt,numZ,dz;
	double emitX,betaX,alphaX,gamX,aveX2,aveXPrime2,aveCrsX,xPrime;
	double emitY,betaY,alphaY,gamY,aveY2,aveYPrime2,aveCrsY,yPrime;
   double x,y,z,px,py,pz,gamma;
   double ax,bx,gx,ay,by,gy,a0x,b0x,g0x,a0y,b0y,g0y;
   double *data;
   FILE *out,*in;
   char fileName[100],fileName1[100],dataName[100],dataSet[100];
   int myrank, nTasks;
   MPI_Status status;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


   recv=(int *)malloc(nTasks*sizeof(int ));

   if(argc < 4)   {
     if(myrank==0) {
       printf("mpirun -np [ ] hdf_track species step maxZ numZ\n");
       //printf("mpirun -np [ ] hdf_track species step maxZ numZ\n");
     } else ;
     MPI_Finalize();  
     exit(0);
   }  else;

   species=atoi(argv[1]);
   step=atoi(argv[2]);
   maxZ=atof(argv[3]);
   numZ=atof(argv[4]);
   dz=maxZ/numZ; 

   sprintf(fileName,"%dParticle%d.h5",species,step);
   if(fopen(fileName,"r")==NULL)  {
     printf("%s is not exited.\n",fileName);
     exit(0);
   } else ;

   sprintf(dataName,"totalCnt");
   if(myrank==0)  
     restoreIntMeta(fileName,dataName,&totalCnt,1);
   else ;
   MPI_Barrier(MPI_COMM_WORLD);
//       MPI_Bcast(&totalCnt,nTasks,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&totalCnt,1,MPI_INT,0,MPI_COMM_WORLD);

   calSub(totalCnt,&cntSub,&start);
   data = (double *)malloc(cntSub*9*sizeof(double ));
   sprintf(dataName,"%d",species);
   restoreData(fileName,dataName,totalCnt,cntSub,start,data,0,9);

	emitX=0.0; betaX=0.0; alphaX=0.0; gamX=0.0;
	emitY=0.0; betaY=0.0; alphaY=0.0; gamY=0.0;
   aveX2=aveXPrime2=aveCrsX=0.0;
   aveY2=aveYPrime2=aveCrsY=0.0;
   aveGam=0.0;
   cnt=0;


   sprintf(fileName1,"%dParticle%d",species,step);
	out=fopen(fileName1,"w");
   for(n=0; n<cntSub; n++)  {
     pz=data[n*dataNum+3];
//     if(pz>minPz && pz<maxPz) {
       z=data[n*dataNum+0];
       x=data[n*dataNum+1];
       y=data[n*dataNum+2];
       pz=data[n*dataNum+3];
       px=data[n*dataNum+4];
       py=data[n*dataNum+5];
	    gamma=sqrt(1+px*px+py*py+pz*pz);
		 fprintf(out,"%g %g %g %g %g %g\n",z,x,y,pz,px,py);

		 xPrime=px/pz;
		 aveX2+=x*x;
		 aveXPrime2+=xPrime*xPrime;
		 aveCrsX+=x*xPrime;

		 yPrime=py/pz;
		 aveY2+=y*y;
		 aveYPrime2+=yPrime*yPrime;
		 aveCrsY+=y*yPrime;

		 aveGam+=gamma;
       cnt+=1.0;
//     }
	}
	fclose(out);
	printf("%s is made.\n",fileName1);

   emitX=sqrt((aveX2*aveXPrime2-aveCrsX*aveCrsX)/cnt/cnt);
   betaX=aveX2/cnt/emitX;
   gamX=aveXPrime2/cnt/emitX;
   alphaX=-aveCrsX/cnt/emitX;

   emitY=sqrt((aveY2*aveYPrime2-aveCrsY*aveCrsY)/cnt/cnt);
   betaY=aveY2/cnt/emitY;
   gamY=aveYPrime2/cnt/emitY;
   alphaY=-aveCrsY/cnt/emitY;

   aveGam/=cnt;
   printf("aveGam=%g, emitX=%g, betaX=%g, gamX=%g, alphaX=%g, emitX_norm=%g\n",aveGam,emitX,betaX,gamX,alphaX,emitX*aveGam);
   printf("aveGam=%g, emitY=%g, betaY=%g, gamY=%g, alphaY=%g, emitY_norm=%g\n",aveGam,emitY,betaY,gamY,alphaY,emitY*aveGam);

   sprintf(fileName1,"twiss");
	out=fopen(fileName1,"w");
   b0x=betaX; a0x=alphaX; g0x=gamX;
   b0y=betaY; a0y=alphaY; g0y=gamY;
   for(z=0; z<maxZ; z+=dz) {
      bx=b0x-2*a0x*z+z*z*g0x;
      ax=a0x-g0x*z;
		gx=g0x;

      by=b0y-2*a0y*z+z*z*g0y;
      ay=a0y-g0y*z;
		gy=g0y;

      fprintf(out,"%g %g %g %g %g %g %g %g %g\n",z,sqrt(emitX*bx),sqrt(emitY*by),bx,ax,gx,by,ay,gy);

      b0x=by; a0x=ax; g0x=gx;
      b0y=by; a0y=ay; g0y=gy;
	}
   fclose(out);
	printf("%s is made.\n",fileName1);

   free(data);
   free(recv);
   MPI_Finalize();
}

void restore1Data(char *fileName,char *dataName,int totalCnt,int cntSub,int start,double *data,int startC)
{
   hid_t file_id,dset_id,plist_id;
   hid_t filespace,memspace;
   hsize_t dimsf[2],count[2],offSet[2],block[2],stride[2];
   herr_t ierr;

   //open file
   plist_id=H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
   file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
   H5Pclose(plist_id);

   //set dataset
   dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
   //file space
   dimsf[0]=totalCnt;
   dimsf[1]=9;
   filespace=H5Screate_simple(2,dimsf,NULL);
   //memory space
   dimsf[0]=cntSub;
   dimsf[1]=1;
   memspace=H5Screate_simple(2,dimsf,NULL);

   stride[0]=1;   stride[1]=1;
   count[0]=1;    count[1]=1;

   //hyperslab in file space
   block[0]=cntSub;  block[1]=1;
   offSet[0]=start;  offSet[1]=startC;
   H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offSet,stride,count,block);
   //hyperslab in memory space
   offSet[0]=0;  offSet[1]=0;
   H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offSet,stride,count,block);
   //read data
   plist_id=H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
   ierr=H5Dread(dset_id,H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);
   H5Dclose(dset_id);
   H5Pclose(plist_id);
   H5Sclose(memspace);
   H5Sclose(filespace);
   H5Fclose(file_id);
}
//lala
void restoreData(char *fileName,char *dataName,int totalCnt,int cntSub,int start,double *data,int startC,int columns)
{
   hid_t file_id,dset_id,plist_id;
   hid_t filespace,memspace;
   hsize_t dimsf[2],count[2],offSet[2],block[2],stride[2];
   herr_t ierr;

   //open file
   plist_id=H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
   file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
   H5Pclose(plist_id);

   //set dataset
   dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);

   //file space
   dimsf[0]=totalCnt;
   dimsf[1]=9;
   filespace=H5Screate_simple(2,dimsf,NULL);
   //memory space
   dimsf[0]=cntSub;
   dimsf[1]=columns;
   memspace=H5Screate_simple(2,dimsf,NULL);

   stride[0]=1;   stride[1]=1;
   count[0]=1;    count[1]=1;

   //hyperslab in file space
   block[0]=cntSub;  block[1]=columns;
   offSet[0]=start;  offSet[1]=startC;
   H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offSet,stride,count,block);
   //hyperslab in memory space
   offSet[0]=0;      offSet[1]=0;
   H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offSet,stride,count,block);
   //read data
   plist_id=H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
   ierr=H5Dread(dset_id,H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);

   H5Pclose(plist_id);
   H5Sclose(memspace);
   H5Sclose(filespace);
   H5Dclose(dset_id);
   H5Fclose(file_id);
}


void calSub(int totalCnt,int *cntSub,int *start)
{
   int i,sub,remain,rank,tmp,*recv,subCnt;;
   int myrank, nTasks;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   recv=(int *)malloc(nTasks*sizeof(int ));

   sub=totalCnt/nTasks;
   remain=totalCnt%nTasks;
   for(rank=0; rank<nTasks; rank++) {
     if(rank<remain)  tmp=sub+1;
     else             tmp=sub;
     if(myrank==rank)  subCnt=tmp; else;
   }
   *cntSub=subCnt;
   MPI_Gather(&subCnt,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);

   tmp=0;
   for(i=0; i<myrank; i++) tmp+=recv[i];
   *start=tmp;
  
   free(recv);
}

void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

double whatMass(int atomNum)
{
  float eMassU;

  eMassU=5.485799e-4;
  if(atomNum==0)        return 1.0;
  else if(atomNum==1)   return 1.0/eMassU;
  else if(atomNum==2)   return 4.0/eMassU;
  else if(atomNum==6)   return 12.0/eMassU;
  else { printf("no list in atomNum data\n"); exit(0); }
}

