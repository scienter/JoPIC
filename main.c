#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include <time.h>

int main(int argc, char *argv[])
{
    int i,j,k,n,s,iteration=0,filter,boost,filterStep,labSaveStep;
    int rnk,suddenDump=OFF;
    double factor,time_spent;
    clock_t begin,end;
    double t;
    char name[100];
    FILE *out;
    Domain D;  
    LaserList *L;
    LoadList *LL;
    External Ext;
    UPML UPml;
    DPML DPml;
    int myrank, nTasks;
    MPI_Status status; 

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    begin=clock();

    if(argc < 2) 
    {  
      printf("mpirun -np N show [inputFile] [dumpNum]\n"); 
      exit(0); 
    }
    if(FindParameters("Domain",1,"filter",argv[1],name)) filter=atoi(name);
    else  filter=0;
    if(FindParameters("Domain",1,"filterStep",argv[1],name)) filterStep=atoi(name);
    else  filterStep=10;

    //parameter setting
    parameterSetting(&D,&Ext,argv[1]);
    if(argc >= 3)
      D.dumpStep = atoi(argv[2]);
    else	;
    MPI_Barrier(MPI_COMM_WORLD);

    //create mesh
    boundary(&D,&Ext);
    MPI_Barrier(MPI_COMM_WORLD);
/*
    //load boost frame's laser
    iif(D.boostOn==1) {
      L=D.laserList;
      while(L->next)  {
        boostLoadLaser2D(&D,L);  
        L=L->next;
      }
      MPI_TransferF_DSX_Yminus(&D,D.numShareDn);
      MPI_TransferF_DSX_Yplus(&D,D.numShareUp);
    }
*/

    //load plasma or load dump file
    if(argc >= 3)
    {   
      iteration=D.dumpStep;
      MPI_Barrier(MPI_COMM_WORLD);
      restoreData(&D,iteration);
//      if(D.saveMode==TXT)
//        restoreData(&D,iteration);
//      else if(D.saveMode==HDF)
//        restoreDumpHDF(&D,iteration);
//      else 	;

      t=D.dt*iteration;
    }
    else 
    {
      LL=D.loadList;
      s=0;
      while(LL->next)
      {
        loadPlasma(&D,LL,s,iteration);
        LL=LL->next;
        s++;
      }
      t=0;
    }

    //rooping time 
    while(iteration<=D.maxStep)
    {

      if(D.tracking==ON && iteration%D.trackSaveStep==0)
      {
        if(D.dimension==2)
          trackID(&D,iteration,D.istart,D.iend,D.jstart,D.jend,0,1);
        else if(D.dimension==3)
          trackID(&D,iteration,D.istart,D.iend,D.jstart,D.jend,D.kstart,D.kend);
      }
      
       //calculating running time     
       end=clock();
       time_spent=(end-begin)/CLOCKS_PER_SEC/60.0;
       if(myrank==0)
       {
         if(time_spent>D.maxTime)
           suddenDump=ON;
         else	;       
       } 
       else	;
       if(myrank==0)
       {
         for(rnk=1; rnk<nTasks; rnk++)
           MPI_Send(&suddenDump,1,MPI_INT,rnk,myrank,MPI_COMM_WORLD);
       }
       else
         MPI_Recv(&suddenDump,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
       MPI_Barrier(MPI_COMM_WORLD);
       if(suddenDump==ON)
       {
         if(D.saveMode==TXT) 	saveDump(D,iteration);
         else if(D.saveMode==HDF) 	saveDumpHDF(D,iteration);
         iteration=D.maxStep+1;
       }
       else	;


       //save File      
       if(iteration%100==0)  {
         saveCenterField(&D,iteration);
         saveCenterDensity(&D,iteration);
       }

       if(iteration%D.saveStep==0 && iteration>=D.saveStart)   
       {
          if(D.fieldSave==ON) { 
            if(D.saveMode==TXT)  	saveField(&D,iteration);
            else if(D.saveMode==HDF)   	saveFieldHDF(D,iteration);
            else 	;
            if(myrank==0)
              printf("field%d is made.\n",iteration);  
          }
          else	;

          if(D.ramanSave==ON) {
	    if(D.saveMode==TXT)		saveRaman(&D,iteration);
            else if(D.saveMode==HDF)	saveRamanHDF(D,iteration);
            else ;
            if(myrank==0)
              printf("raman%d is made.\n",iteration);  
          }
          else	;

          if(D.particleSave==ON) { 
             saveParticle(&D,iteration);
//            if(D.saveMode==TXT)		saveParticle(&D,iteration);
//            else if(D.saveMode==HDF)	saveParticleHDF(D,iteration,D.nSpecies);
//            else ;
            if(myrank==0)
              printf("particle%d is made.\n",iteration); 
          }
          else	;

          if(D.rhoSave==ON) { 
            if(D.saveMode==TXT)		saveRho(&D,iteration);
            else if(D.saveMode==HDF)	saveDensityHDF(D,iteration,D.nSpecies);
            else			;
            if(myrank==0)
              printf("rho%d is made.\n",iteration);    
          }
          else	;


          if(D.dumpSave==ON && iteration>=D.dumpStart) { 
              saveDump(D,iteration);
//            if(D.saveMode==TXT) 	saveDump(D,iteration);
//            else if(D.saveMode==HDF) 	saveDumpHDF(D,iteration);
//            else 	;
            if(myrank==0)
              printf("dump%d is made.\n",iteration);              
          }

          if(D.probeNum>0) { 
//            saveProbe(&D,iteration);
//            if(myrank==0)
//              printf("probe%d is made.\n",iteration);              
          }

       }



//       if(nTasks==1)  periodY1core(&D);   
//       else           periodY(&D);

       fieldSolveC(&D);
       fieldTransferC(&D);
       fieldSolve(&D);

       //load laser
       if(D.boostOn==OFF)
       {
         L=D.laserList;
         while(L->next)  {
           loadLaser(&D,L,t); 
//         if(L->direction==1)     loadLaser2D(&D,L,t); 
//         else if(L->direction==-1)     loadLaserOpp2D(&D,L,t); 
           L=L->next;
         }
       }
       fieldTransfer(&D);

       interpolation(&D,&Ext);
       particlePush(&D);

       updateCurrent(&D);

       if (iteration>=D.nx && D.moving==ON && D.boostOn==OFF)
       {
         movingDomain(&D);
         LL=D.loadList;
         s=0;
         while(LL->next)
         {
           loadMovingPlasma(&D,LL,s,iteration);
           LL=LL->next;
           s++;
         }
         rearrangeParticles(&D);
         if(D.M>1)
           particleShareY(&D);
         if(D.N>1)
           particleShareZ(&D);
         removeEdge(&D);
       }
       else if(D.boostOn==ON) {
         LL=D.loadList;
         s=0;
         while(LL->next)
         {
           loadMovingPlasma(&D,LL,s,iteration);
           LL=LL->next;
           s++;
         }
         movingDomain(&D);
         rearrangeParticles(&D);
         if(D.M>1)
           particleShareY(&D);
         removeEdge(&D);
       }
       else
       {
          rearrangeParticles(&D);
          if(D.M>1)
            particleShareY(&D);
          if(D.N>1)
            particleShareZ(&D);
          removeEdge(&D);
       }

       //time update
       if(iteration%10==0 && myrank==0)  
          printf("iteration = %d\n",iteration);           
       iteration+=1;
       t=D.dt*iteration;  

    }     //end of time roop                  

    if(D.tracking==ON)
      saveTracking(&D);

    end=clock();
    time_spent=(end-begin)/CLOCKS_PER_SEC;

    //make 'report' file
    if(myrank==0)
    {
      sprintf(name,"report");
      out = fopen(name,"w");
      fprintf(out,"nx=%d\n",D.nx);
      fprintf(out,"ny=%d\n",D.ny);
      fprintf(out,"nz=%d\n",D.nz);
      fprintf(out,"cores=%d\n",nTasks);
      fprintf(out,"running time=%gm\n",time_spent/60.0);
      fclose(out);
    }
    else	;

    cleanMemory(&D);
    
    MPI_Finalize();

    return 0;
}
