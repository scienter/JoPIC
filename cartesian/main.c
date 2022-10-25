#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>
#include <time.h>

void saveDefParticle(Domain *D,int iteration);

int main(int argc, char *argv[])
{
  	int i,j,k,n,s,iteration=0,filter,boost,filterStep,labSaveStep;
  	int rnk,suddenDump=OFF,startI,endI,half,non_half;
  	double x,factor,time_spent,t;
  	clock_t begin,end;
  	struct tm *t_now;
  	time_t timer; 	//measure time
  	char name[100];
  	FILE *out;
  	Domain D;  
  	LaserList *L;
  	LoadList *LL;
  	PlasmaLens *PL;	 
  	External Ext;
  	UPML UPml;
  	DPML DPml;
  	int myrank, nTasks,rankX;
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

  	timer=time(NULL);
  	t_now=localtime(&timer);
  	sprintf(name,"report");
  	if(myrank==0) {  
  	  	out = fopen(name,"a");
  	  	fprintf(out,"simulation start.\n");
  	  	fprintf(out,"%d-%d-%d %d:%d:%d\n",t_now->tm_year+1900,t_now->tm_mon+1,t_now->tm_mday,t_now->tm_hour,t_now->tm_min,t_now->tm_sec);
  	} else ;

  	//parameter setting
  	parameterSetting(&D,&Ext,argv[1]);
  	if(argc >= 3)      D.dumpStep = atoi(argv[2]); else;
  	if(myrank==0)     saveDensityProfile(&D); else ;
	MPI_Barrier(MPI_COMM_WORLD);

  	//create mesh
  	boundary(&D,&Ext);
  	MPI_Barrier(MPI_COMM_WORLD);

  	//load plasma or load dump file
  	if(argc >= 3)  {   
   	iteration=D.dumpStep;
		//     assignDefParticle(&D);
    	restoreDump(&D,iteration);
    	t=D.dt*iteration;

    	LL=D.loadList;
    	while(LL->next)      {
      	LL->createGuard=D.maxX+LL->RDef*2.0;
      	LL->deleteGuard=D.maxX-LL->RDef*2.0;        
      	LL=LL->next;
   	}
		//      deleteDefParticle(&D);
	}  else   {
   	assignDefParticle(&D);
   	LL=D.loadList;
   	s=0;
    	while(LL->next)      {
    	  	loadPlasma(&D,LL,s,D.istart,D.iend,iteration);
    	  	if(LL->pair==ON) solveCharge(&D,LL,D.RhoPair,s,D.istart,D.iend,-1.0,half=0,non_half=1); else ;
	 	  	// loading track particle for loading Plasma
    		if(D.tracking==ON) track(&D,D.istart,D.iend,s); else ;

    	  	LL->createGuard=D.maxX+LL->RDef*2.0;
    	  	LL->deleteGuard=D.maxX-LL->RDef*2.0;        
    	  	LL=LL->next;
    	  	s++;
    	}
    	deleteDefParticle(&D);
    	t=0;

    	//pair charge
    	D.shareF[0]=D.RhoPair;
    	MPI_TransferJ_Xplus(&D,1);
    	MPI_TransferJ_Xminus(&D,1);
    	MPI_TransferJ_Yplus(&D,1);
    	MPI_TransferJ_Yminus(&D,1);
    	MPI_TransferJ_Zplus(&D,1);
    	MPI_TransferJ_Zminus(&D,1);
    	D.RhoPair=D.shareF[0];
  	}
  	MPI_Barrier(MPI_COMM_WORLD);



  	// load Laser Shot
  	L=D.laserList;
  	while(L->next)  {
  	  	if(L->loadMethod==Shot) shotLaser(&D,L); else ;
  	  	L=L->next;
  	}


  	//rooping time 
  	while(iteration<=D.maxStep)
  	{

 		solveF(&D,D.Pr,D.Pl,D.Sr,D.Sl,D.EzNow,D.BzNow,half=0,non_half=1);
      //if(myrank==0) printf("iteration=%d, solveF is done.\n",iteration);
		 		  
		//     //calculating running time     
		//       end=clock();
		//       time_spent=(end-begin)/CLOCKS_PER_SEC/60.0;
		//       if(myrank==0)
		//       {
		//         if(time_spent>D.maxTime)
		//           suddenDump=ON;
		//         else	;       
		//       } 
		//       else	;
		//       if(myrank==0) {
		//         for(rnk=1; rnk<nTasks; rnk++)
		//           MPI_Send(&suddenDump,1,MPI_INT,rnk,myrank,MPI_COMM_WORLD);
		//       } else;
		//       MPI_Barrier(MPI_COMM_WORLD);
		//       if(suddenDump==ON) {
		//         saveDump(D,iteration);
		//         iteration=D.maxStep+1;
		//       } else	;

    	// Beam loading
		LL=D.loadList;
		if(iteration==LL->loadingStep) {      
      	s=0;
     		while(LL->next)      {
      		loadBeam(&D,LL,s,iteration);
				if(D.tracking==ON) track(&D,D.istart,D.iend,s); else ;
      	  	LL=LL->next;
      	  	s++;
     		}
			if(myrank==0) printf("At step%d, beam loading is done.\n",iteration); else ;
		} else ;
		MPI_Barrier(MPI_COMM_WORLD);
      //if(myrank==0) printf("iteration=%d, beam loading is done.\n",iteration);


   	//resolChange
   	if(D.resolChange==ON)  {
      	if(iteration==D.resolStep)   {
        		saveBDump(D,iteration);
        		saveEDump(D,iteration);
        		saveJDump(D,iteration);
				//           if(D.resolHigh==ON)
				//             saveDumpParticleResolHDF(&D,iteration);
				//           else ;
      	} else	;
      	if(iteration-1==D.resolStep)   {
        		saveBDump(D,iteration);
        		saveJDump(D,iteration);
				//           if(D.resolLow==ON)
				//             saveDumpParticleResolHDF(&D,iteration);
				//           else ;
      	} else	;
      	if(iteration+1==D.resolStep)   {
        		saveBDump(D,iteration);
        		saveJDump(D,iteration);
      	} else	;
    	}      else; 

   	//save File      
    	if(iteration%D.saveStep==0 && iteration>=D.saveStart)   {
      	saveFile(D,iteration);
      	end=clock();
      	time_spent=(end-begin)/CLOCKS_PER_SEC/60.0;
      	if(myrank==0) {
        		fprintf(out,"Time duration at %dth saveFile:%.4gmin\n",iteration,time_spent);
        		printf("Time duration at %dth saveFile:%.4gmin\n",iteration,time_spent);
      	} else ;
    	}  else	;
    	//save dump File
    	if(D.dumpSave==ON && iteration>=D.dumpStart && iteration%D.dumpSaveStep==0) 
      	saveDump(D,iteration); 
    	else	;

   	//save center field      
   	if(iteration%D.centerStep==0)  {
      	saveCenterField(&D,iteration);
			//         saveCenterDensity(&D,iteration);
    	}  else	;
      //if(myrank==0) printf("iteration=%d, save Center is done.\n",iteration);

    	if(D.consCheck==ON) checkEnergyCons(&D,iteration); else ;
      //if(myrank==0) printf("iteration=%d, check energy conservation.\n",iteration);

    	// tracking
		particleTracking(&D,iteration);
    	MPI_Barrier(MPI_COMM_WORLD);
    	if(D.tracking==ON && (iteration+1)%100==0) saveTrack(&D,iteration); else ;
      //if(myrank==0) printf("iteration=%d, particle tracking is done.\n",iteration);


    	interpolation(D,&Ext);
		//if(myrank==0) printf("iteration=%d, interpolation\n",iteration);

    	// Plasma lens
    	PL=D.lensList;
    	while(PL->next)      {
      	plasmaLens(&D,PL,iteration);
      	PL=PL->next;
    	}
		//if(myrank==0) printf("iteration=%d, plasma lens is done\n",iteration);

    	particlePush(&D,iteration);
		//if(myrank==0) printf("iteration=%d, push\n",iteration);

    	if(D.fieldIonizationONOFF==ON) fieldIonization(D,iteration); else;
    	MPI_Barrier(MPI_COMM_WORLD);
		//if(myrank==0) printf("iteration=%d, fieldInization\n",iteration);

    	if(D.fieldType==Split)	updateCurrent_Split(D,iteration);
	  	else 							updateCurrent(D);
		//if(myrank==0) printf("iteration=%d, current\n",iteration);

		// Moving Domain
   	x=D.movingV*iteration+D.shiftStart;	// 'x' is a guard for moving domain.
   	if(D.moving==ON && x>D.maxXDomain-1)    
		{

      	movingDomain(&D,iteration);
			//if(myrank==0) printf("iteration=%d, moving domain\n",iteration);

      	if(myrank/(D.M*D.N)==D.L-1) {
		   	LL=D.loadList;
        		s=0;
        		while(LL->next)      {
          		loadPlasma(&D,LL,s,D.iend-1,D.iend,iteration);
			   	if(D.tracking==ON) track(&D,D.iend-1,D.iend,s); else ;
			   	LL=LL->next;
          		s++;
        		}
      	} else ;

			//	 assignDefParticle(&D);
			//         if(iteration%D.saveStep==0 && iteration>=D.saveStart)   {
			//           if(myrank==0) saveDefParticle(&D,iteration); else ;
			//        } else ;

      	deleteDefParticle(&D);
			//if(myrank==0) printf("iteration=%d, deleteDefParticle\n",iteration);

      	rearrangeParticles(&D,iteration);
      	particleShareX(D);
      	particleShareY(D);
      	particleShareZ(D);

      	if(myrank/(D.M*D.N)==D.L-1) {
		   	LL=D.loadList;
        		s=0;
        		while(LL->next)      {
         		if(LL->pair==ON) solveCharge(&D,LL,D.RhoPair,s,D.iend-1,D.iend,-1.0,half=0,non_half=1); else ;
			   	LL=LL->next;
          		s++;
        		}
      	} else ;
      	D.shareF[0]=D.RhoPair;
      	MPI_TransferJ_Yplus(&D,1);
      	MPI_TransferJ_Yminus(&D,1);
      	MPI_TransferJ_Zplus(&D,1);
      	MPI_TransferJ_Zminus(&D,1);
      	D.RhoPair=D.shareF[0];

      	removeEdge(D);
			//if(myrank==0) printf("iteration=%d, removeEdge\n",iteration);

   	} 
	
		// Non Moving Domain
		else       
		{
      	rearrangeParticles(&D,iteration);
      	particleShareX(D);
      	particleShareY(D);
      	particleShareZ(D);
			//if(myrank==0) printf("iteration=%d, rearrange\n",iteration);

      	removeEdge(D);
			//printf("nomoving, iteration=%d, removeEdge\n",iteration);
    	}
          

		fieldSolve1(D,t,iteration);
		//if(myrank==0 && iteration>=500) printf("iteration=%d, fieldsolve1\n",iteration);

 		solveF(&D,D.PrC,D.PlC,D.SrC,D.SlC,D.Ez,D.Bz,half=1,non_half=0);

   	fieldSolve2(D,t,iteration);
		//if(myrank==0) printf("iteration=%d, fieldSolve2 is done.\n",iteration);       

    	//time update
    	if(iteration%10==0 && myrank==0) printf("iteration = %d\n",iteration);
    	else ;
    	iteration+=1;
    	t=D.dt*iteration;  
  	}     //end of time roop                  


	end=clock();
  	time_spent=(end-begin)/CLOCKS_PER_SEC;

  	//make 'report' file
  	if(myrank==0) {
   	fprintf(out,"nx=%d, ",D.nx);
   	fprintf(out,"ny=%d, ",D.ny);
   	fprintf(out,"nz=%d\n",D.nz);
   	fprintf(out,"cores=%d, \n",nTasks);
   	fprintf(out,"nSpecies=%d\n",D.nSpecies);
   	fprintf(out,"running time=%.4gm\n",time_spent/60.0);
   	fprintf(out,"\n");
   	fclose(out);
  	}  else ;

  	cleanMemory(&D);

  	MPI_Finalize();

  	return 0;
}
