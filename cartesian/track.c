#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>


void saveDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void share(Domain *D,double *data,int dataNum);
void saveFieldComp_serial(Domain *D,char *fileName,int iteration);


void saveTrack(Domain *D,int iteration)
{
	int dataNum,n,i,rank,minXSub,minYSub,minZSub,istart,jstart,kstart,MAX;
	double dt,unitX,unitY,unitZ,*recv;
	char fileName[100];

	int nTasks,myrank;
   MPI_Status status;	
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   hid_t file_id,group_id,dset_id,filespace;
   herr_t hdf_status;

	istart=D->istart;	jstart=D->jstart;	kstart=D->kstart;
	minXSub=D->minXSub; minYSub=D->minYSub; minZSub=D->minZSub;
	MAX=D->trackStep;

	dataNum=MAX*D->idNums*6;
	share(D,D->trackParticle,dataNum);

	sprintf(fileName,"track.h5");
	if(myrank==0) {
		saveFieldComp_serial(D,fileName,iteration);
//      saveFieldComp_serial(D->trackX,fileName,"trackX",D->idNums,D->maxStep,MAX,iteration);
//      saveFieldComp_serial(D->trackY,fileName,"trackY",D->idNums,D->maxStep,MAX,iteration);
//      saveFieldComp_serial(D->trackZ,fileName,"trackZ",D->idNums,D->maxStep,MAX,iteration);
//      saveFieldComp_serial(D->trackUx,fileName,"trackUx",D->idNums,D->maxStep,MAX,iteration);
//      saveFieldComp_serial(D->trackUy,fileName,"trackUy",D->idNums,D->maxStep,MAX,iteration);
//      saveFieldComp_serial(D->trackUz,fileName,"trackUz",D->idNums,D->maxStep,MAX,iteration);

		printf("Intermediate tracking is done at %d\n",iteration);
	} else ;

	for(i=0; i<dataNum; i++) {
		D->trackParticle[i]=0.0;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	

}


void particleTracking(Domain *D,int iteration)
{
	int i,j,k,istart,iend,jstart,jend,kstart,kend;
	int index,minXSub,minYSub,minZSub,idNums,tIdx,MAX;
	double unitX,unitY,unitZ;
	FILE *out;
	char fileName[100];
	ptclList *p;

	istart=D->istart;    iend=D->iend;
	jstart=D->jstart;    jend=D->jend;
	kstart=D->kstart;    kend=D->kend;
	minXSub=D->minXSub; minYSub=D->minYSub; minZSub=D->minZSub;
	idNums=D->idNums;
	MAX=D->trackStep;

	unitX=D->dx*D->lambda;
	unitY=D->dy*D->lambda;
	unitZ=D->dz*D->lambda;
	tIdx=iteration%MAX;
	for(i=istart; i<iend; i++)
		for(j=jstart; j<jend; j++)
			for(k=kstart; k<kend; k++)   {
				p=D->track[i][j][k].head[0]->pt;
				while(p) {
					index=(int)(p->p1Old1);
					D->trackParticle[index*(MAX*6)+tIdx*6+0]=(0.5*(p->x+i+p->oldX)-istart+minXSub)*unitX;
					D->trackParticle[index*(MAX*6)+tIdx*6+1]=(0.5*(p->y+j+p->oldY)-jstart+minYSub)*unitY;
					D->trackParticle[index*(MAX*6)+tIdx*6+2]=(0.5*(p->z+k+p->oldZ)-kstart+minZSub)*unitZ;
					D->trackParticle[index*(MAX*6)+tIdx*6+3]=p->p1;
					D->trackParticle[index*(MAX*6)+tIdx*6+4]=p->p2;
					D->trackParticle[index*(MAX*6)+tIdx*6+5]=p->p3;
/*
					sprintf(fileName,"%d",index);
					out=fopen(fileName,"a");
					fprintf(out,"%.20g %.20g %g %.20g %.20g %g\n",
						(p->x+i-istart+minXSub)*unitX,
						(p->y+j-jstart+minYSub)*unitY,
						(p->z+k-kstart+minZSub)*unitZ,
						p->p1,
						p->p2,
						p->p3);
					fclose(out);
*/
					p=p->next;
				}
			}
 
}

void track(Domain *D,int startI,int endI,int s)
{
	int i,j,k,istart,iend,jstart,jend,kstart,kend,cnt,FLAG,tmpInt;
	double id,core;
    Particle ***particle;
    particle=D->particle;
	IDCore *tmpTr,*prevTr,*idcore;
	ptclList *New,*p;
    int myrank, nTasks;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    for(i=startI; i<endI; i++)  
	{
    	for(j=jstart; j<jend; j++)	
		{
        	for(k=kstart; k<kend; k++)   
		  	{
            	p=particle[i][j][k].head[s]->pt;     
            	while(p)    {    
					id=p->index; core=p->core;

					cnt=1;
					idcore=D->trhead->pt;
					while(idcore) {
						if(cnt==1) prevTr=idcore; else ;

						FLAG=0;
						if(id==idcore->id) {
                     New = (ptclList *)malloc(sizeof(ptclList));
				        	New->next = D->track[i][j][k].head[0]->pt;
							D->track[i][j][k].head[0]->pt = New;	
							
							New->x=p->x; New->oldX=p->oldX;
							New->y=p->y; New->oldY=p->oldY;
							New->z=p->z; New->oldZ=p->oldZ;
							New->E1=p->E1; New->E2=p->E2; New->E3=p->E3;
							New->B1=p->B1; New->B2=p->B2; New->B3=p->B3;
							New->p1=p->p1; New->p2=p->p2; New->p3=p->p3;
							New->p1Old1=p->p1Old1; New->p2Old1=p->p2Old1; New->p3Old1=p->p3Old1;
							New->p1Old1=idcore->index; 
							New->weight=p->weight; New->charge=p->charge;
							New->index=p->index;	New->core=p->core;
							FLAG=1;
						} else ;

						if(FLAG==1) {
							if(cnt==1) {
								tmpTr=idcore->next;
								D->trhead->pt=tmpTr;
								idcore->next=NULL;
								free(idcore);
								idcore=D->trhead->pt;
								cnt=1;
							} else {
								prevTr->next = idcore->next;
								idcore->next=NULL;
								free(idcore);
								idcore=prevTr->next;
								cnt++;
							}
							break;
						} else {
							prevTr=idcore;
							idcore=idcore->next;
							cnt++;
						}

					}
            		p=p->next;
          		}		//End of while(p)
				  
        	}	//End of for(k)
		}		//End of for(j)
	 }			//End of for(i)
}

void saveFieldComp_serial(Domain *D,char *fileName,int iteration)
{
    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,memspace,dataspace,ticspace;
    hsize_t dimsf[2],dimsm[2],count[2],offset[2],block[2],stride[2];
	 double *data;
	 int MAX,i,tIdx;
	 char dataName[100];

    data=(double *)malloc((D->trackStep*6)*sizeof(double ));
    file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
	
	MAX=D->trackStep;
	tIdx=(iteration-(MAX-1))%MAX;

	dimsf[0]=D->maxStep;
   dimsf[1]=6;
   dataspace=H5Screate_simple(2,dimsf,NULL);

	dimsm[0]=MAX;
	dimsm[1]=6;
	memspace=H5Screate_simple(2,dimsm,NULL);
	
	count[0]=MAX;
	count[1]=6;

	stride[0]=1;
	stride[1]=1;

	block[0]=1;
	block[1]=1;

	for(i=0; i<D->idNums; i++) {
		offset[0]=iteration-(MAX-1);
		offset[1]=0;

		memcpy(data,&D->trackParticle[i*(MAX*6)+tIdx*6],MAX*6*sizeof(double));
		sprintf(dataName,"%d",i);
		dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
//    dataspace = H5Dget_space(dset_id);
   	status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, stride, count, block);
   	status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,dataspace,H5P_DEFAULT,data);
	}

	 H5Dclose(dset_id);
	 H5Sclose(memspace);
	 H5Sclose(dataspace);
	 H5Fclose(file_id);
}

void share(Domain *D,double *data,int dataNum) 
{
	int i,rank;
	double *recv;

	int myrank, nTasks;
   MPI_Status status;	
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   recv=(double *)malloc(dataNum*sizeof(double ));

   if(myrank!=0)
      MPI_Send(data,dataNum,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
   else {
		for(rank=1; rank<nTasks; rank++) {
			MPI_Recv(recv,dataNum,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
			for(i=0; i<dataNum; i++)   data[i]+=recv[i];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
//   MPI_Bcast(data,dataNum,MPI_DOUBLE,0,MPI_COMM_WORLD);    	
	free(recv);
}
