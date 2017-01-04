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
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz,bias,biasY,biasZ,offSetY;
    int nxSub,nySub,nzSub,cnt,totalCnt,index,start,ii,tmp,rankY,rankZ;
    int offset[3];
    char name[100],name2[100];
    double *saveDouble;
    int *saveInt;
    void saveIntArray();
    Particle ***particle;
    particle=D.particle;
    ptclList *p;
    LoadList *LL;

    int myrank, nTasks;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    int *recv,*offSetRank;  
    recv = (int *)malloc(nTasks*sizeof(int ));
    offSetRank = (int *)malloc(nTasks*sizeof(int ));

    hid_t file_id,group_id,dset_id,filespace;
    hsize_t metaDim[1];
    herr_t status;

    istart=D.istart;
    iend=D.iend;
    jstart=D.jstart;
    jend=D.jend;
    kstart=D.kstart;
    kend=D.kend;
    nxSub=D.nxSub;
    nySub=D.nySub;
    nzSub=D.nzSub;
   
    metaDim[0]=1;
  
    sprintf(name,"dump%d.h5",iteration);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    switch((D.fieldType-1)*3+D.dimension) {
    //2D
    case (Split-1)*3+2:
      ny=D.ny+5;      
      nx=D.nx+5;     
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
      if(myrank==0)
        saveIntMeta(name,"/minXSub",&(D.minXSub));
      else  	;
     
      offset[0]=0;
      offset[1]=D.minYSub-D.minYDomain+bias;
      saveFieldComp_2D(D.Ex,name,"/Ex",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.Pr,name,"/Pr",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.Pl,name,"/Pl",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.Bx,name,"/Bx",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.Sr,name,"/Sr",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.Sl,name,"/Sl",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.ExC,name,"/ExC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.PrC,name,"/PrC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.PlC,name,"/PlC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.BxC,name,"/BxC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.SrC,name,"/SrC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.SlC,name,"/SlC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.Jx,name,"/Jx",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.Jy,name,"/Jy",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.Jz,name,"/Jz",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.JxOld,name,"/JxOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.JyOld,name,"/JyOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_2D(D.JzOld,name,"/JzOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,offset);
      MPI_Barrier(MPI_COMM_WORLD);

      istart=D.istart;
      iend=D.iend;
      jstart=D.jstart;
      jend=D.jend;
      for(s=0; s<D.nSpecies; s++)
      {
        if(myrank==0)
        {  
          file_id=H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);
          sprintf(name2,"%dParticle",s);          
          group_id=H5Gcreate2(file_id,name2,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          H5Gclose(group_id);
          H5Fclose(file_id);
        }
        else 	;
        
        k=0;
        cnt=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              cnt++;
              p=p->next;
            }
          }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(&cnt,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);

        start=0;
        for(i=0; i<myrank; i++)
          start+=recv[i];
        totalCnt=0;
        for(i=0; i<nTasks; i++)
          totalCnt+=recv[i];
        saveDouble = (double *)malloc(cnt*sizeof(double ));
        saveInt = (int *)malloc(cnt*sizeof(int ));

        for(i=0; i<nTasks; i++)
        {
          tmp=0;
          for(ii=0; ii<i; ii++)
            tmp+=recv[ii];
          offSetRank[i]=tmp;
        }              
 
        sprintf(name2,"%dParticle/offSet",s);          
        if(myrank==0)
          saveIntArray(name,name2,offSetRank,nTasks);
        else	;
        sprintf(name2,"%dParticle/counts",s);          
        if(myrank==0)
          saveIntArray(name,name2,recv,nTasks);
        else	;
 
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              saveDouble[index]=p->x;
              index++;
              p=p->next;
            }
          }
        sprintf(name2,"%dParticle/x",s);          
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              saveDouble[index]=p->y;
              index++;
              p=p->next;
            }
          }
        sprintf(name2,"%dParticle/y",s);          
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              saveDouble[index]=p->p1;
              p=p->next;
              index++;
            }
          }
        sprintf(name2,"%dParticle/px",s);          
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              saveDouble[index]=p->p2;
              index++;
              p=p->next;
            }
          }
        sprintf(name2,"%dParticle/py",s);          
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              saveDouble[index]=p->p3;
              index++;
              p=p->next;
            }
          }
        sprintf(name2,"%dParticle/pz",s);          
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              saveInt[index]=p->index;
              index++;
              p=p->next;
            }
          }
        sprintf(name2,"%dParticle/index",s);          
        saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              saveInt[index]=p->core;
              p=p->next;
              index++;
            }
          }
        sprintf(name2,"%dParticle/core",s);          
        saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              saveInt[index]=i;
              p=p->next;
              index++;
            }
          }
        sprintf(name2,"%dParticle/i",s);          
        saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              saveInt[index]=j;
              p=p->next;
              index++;
            }
          }
        sprintf(name2,"%dParticle/j",s);          
        saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
       
        if(myrank==0)
        {
          sprintf(name2,"/%dParticle/totalCnt",s);
          saveIntMeta(name,name2,&(totalCnt));
        }
        else  	;

        free(saveDouble);
        free(saveInt);
      }

      if(myrank==0)
        saveIntMeta(name,"/nSpecies",&(D.nSpecies));
      else  	;           
    break;

    //3D
    case (Split-1)*3+3:
      nz=D.nz+5;      
      ny=D.ny+5;      
      nx=D.nx+5;     
      rankY=myrank%D.M;
      rankZ=myrank/D.M;
      if(rankY==0)  {      
        jstart=0;
        nySub+=2;    
        biasY=0;
      }
      else if(rankY==D.M-1)  {
        jend+=3;
        nySub+=3; 
        biasY=2;
      }
      else   {
        biasY=2; 
      } 
      if(rankZ==0)  {      
        kstart=0;
        nzSub+=2;    
        biasZ=0;
      }
      else if(rankZ==D.N-1)  {
        kend+=3;
        nzSub+=3; 
        biasZ=2;
      }
      else   {
        biasZ=2; 
      } 

      istart=0;
      iend+=3; 
      nxSub+=5;
      if(myrank==0)
        saveIntMeta(name,"/minXSub",&(D.minXSub));
      else  	;
     
      offset[1]=0;
      offset[0]=D.minYSub-D.minYDomain+biasY;
      offset[2]=D.minZSub-D.minZDomain+biasZ;
      saveFieldComp_3D(D.Ex,name,"/Ex",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.Pr,name,"/Pr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.Pl,name,"/Pl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.Bx,name,"/Bx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.Sr,name,"/Sr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.Sl,name,"/Sl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.ExC,name,"/ExC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.PrC,name,"/PrC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.PlC,name,"/PlC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.BxC,name,"/BxC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.SrC,name,"/SrC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.SlC,name,"/SlC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.Jx,name,"/Jx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.Jy,name,"/Jy",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.Jz,name,"/Jz",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.JxOld,name,"/JxOld",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.JyOld,name,"/JyOld",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp_3D(D.JzOld,name,"/JzOld",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);

      istart=D.istart;
      iend=D.iend;
      jstart=D.jstart;
      jend=D.jend;
      kstart=D.kstart;
      kend=D.kend;
      for(s=0; s<D.nSpecies; s++)
      {
        if(myrank==0)
        {  
          file_id=H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);
          sprintf(name2,"%dParticle",s);          
          group_id=H5Gcreate2(file_id,name2,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          H5Gclose(group_id);
          H5Fclose(file_id);
        }
        else 	;
        
        cnt=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++) 
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                cnt++;
                p=p->next;
              }
            }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(&cnt,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);

        start=0;
        for(i=0; i<myrank; i++)
          start+=recv[i];
        totalCnt=0;
        for(i=0; i<nTasks; i++)
          totalCnt+=recv[i];
        for(i=0; i<nTasks; i++)
        {
          tmp=0;
          for(ii=0; ii<i; ii++)
            tmp+=recv[ii];
          offSetRank[i]=tmp;
        }               
        sprintf(name2,"%dParticle/offSet",s);          
        if(myrank==0)
          saveIntArray(name,name2,offSetRank,nTasks);
        else	;
        sprintf(name2,"%dParticle/counts",s);          
        if(myrank==0)
          saveIntArray(name,name2,recv,nTasks);
        else	;

        saveDouble = (double *)malloc(cnt*sizeof(double ));
        saveInt = (int *)malloc(cnt*sizeof(int ));

        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                saveDouble[index]=p->x;
                p=p->next;
                index++;
              }
            }
        sprintf(name2,"%dParticle/x",s);          
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                saveDouble[index]=p->y;
                p=p->next;
                index++;
              }
            }
        sprintf(name2,"%dParticle/y",s);          
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                saveDouble[index]=p->z;
                p=p->next;
                index++;
              }
            }
        sprintf(name2,"%dParticle/z",s);          
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                saveDouble[index]=p->p1;
                p=p->next;
                index++;
              }
            }
        sprintf(name2,"%dParticle/px",s);          
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                saveDouble[index]=p->p2;
                p=p->next;
                index++;
              }
            }
        sprintf(name2,"%dParticle/py",s);          
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                saveDouble[index]=p->p3;
                p=p->next;
                index++;
              }  
            }
        sprintf(name2,"%dParticle/pz",s);          
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                saveInt[index]=p->index;
                p=p->next;
                index++;
              }
            }
        sprintf(name2,"%dParticle/index",s);          
        saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                saveInt[index]=p->core;
                p=p->next;
                index++;
              }
            }
        sprintf(name2,"%dParticle/core",s);          
        saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                saveInt[index]=i;
                p=p->next;
                index++;
              }
            }
        sprintf(name2,"%dParticle/i",s);          
        saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                saveInt[index]=j;
                p=p->next;
                index++;
              }
            }
        sprintf(name2,"%dParticle/j",s);          
        saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                saveInt[index]=k;
                p=p->next;
                index++;
              }
            }
        sprintf(name2,"%dParticle/k",s);          
        saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);
        MPI_Barrier(MPI_COMM_WORLD);
       
        if(myrank==0)
        {
          sprintf(name2,"/%dParticle/totalCnt",s);
          saveIntMeta(name,name2,&(totalCnt));
        }
        else  	;

        free(saveDouble);
        free(saveInt);
      }

      if(myrank==0)
        saveIntMeta(name,"/nSpecies",&(D.nSpecies));
      else  	;           
    break;

  }		//End of switch(dimension....)

  free(recv);
  free(offSetRank);
}

void saveIntArray(char *fileName,char *dataName,int *data,int numData)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[2];
  herr_t status;

  metaDim[0]=numData;
  metaDim[1]=1;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(2,metaDim,NULL);
  dset_id=H5Dcreate(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

