#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void MPI_TransferFNew_Xminus(Domain *D,int numField,int ny,int share)
{
	int i,j,m,num,start,end,numMode,n;
	int istart,iend,jstart,jend;
	int myrank, nTasks,rank; 
	double *data;

	MPI_Status status;        

	istart=D->istart;    iend=D->iend;
	jstart=D->jstart;    jend=D->jend;
	numMode=D->numMode;

	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);           
	
	rank=myrank;

	// Left to Right
	num=numMode*numField*ny*share;
	data = (double *)malloc(num*sizeof(double ));

   //Transferring even ~ odd cores 
	if(rank%2==0 && rank!=D->L-1)
	{
		MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
		start=0; 
		for(n=0; n<numField; n++)
   		for(m=0; m<numMode; m++)		
				for(i=0; i<share; i++) {
	      		for(j=0; j<ny; j++)  
	        			D->shareF[n][m][iend+i][j]=data[start+j]; 
	        		start+=ny; 
				}
	}
	else if(rank%2==1) 
	{
		start=0; 
		for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)		
				for(i=0; i<share; i++) {
		   		for(j=0; j<ny; j++) 
		       		data[start+j]=D->shareF[n][m][i+istart][j]; 
		     		start+=ny; 
		   	}
		MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD); 
	}	
	MPI_Barrier(MPI_COMM_WORLD);

	//Transferring odd ~ even cores                    
	if(rank%2==1 && rank!=D->L-1)
	{
		MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
		start=0; 
		for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)	
				for(i=0; i<share; i++) {
	   			for(j=0; j<ny; j++)
	         		D->shareF[n][m][iend+i][j]=data[start+j]; 
	        		start+=ny; 
	   		}
	}
	else if(rank%2==0 && rank!=0)
	{
		start=0; 
	  	for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)	
	   		for(i=0; i<share; i++) {
	     			for(j=0; j<ny; j++) 
	         		data[start+j]=D->shareF[n][m][i+istart][j]; 
	       		start+=ny; 
	     		}
		MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
	}             
	MPI_Barrier(MPI_COMM_WORLD);
	free(data);
}

void MPI_TransferFNew_Xplus(Domain *D,int numField,int ny,int share)
{
   int i,j,m,num,numMode,n;
   int istart,iend,jstart,jend;
   int myrank, nTasks,rank,start; 
   double *data;

   MPI_Status status;         
   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   numMode=D->numMode;
  
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

   rank=myrank/D->M;

   num=numMode*numField*ny*(share-1);
   data = (double *)malloc(num*sizeof(double ));

  	//Transferring even ~ odd cores 
  	if(rank%2==1)
  	{
  		MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
  		start=0;
  		for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
  		  		for(i=1; i<share; i++) {
  		   		for(j=0; j<ny; j++)
  		      		D->shareF[n][m][istart-i][j]=data[start+j]; 
  		      	start+=ny; 
  		    	}
  	}
  	else if(rank%2==0 && rank!=D->L-1)
  	{
  		start=0; 
  		for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
  		  		for(i=1; i<share; i++) {
  		    		for(j=0; j<ny; j++)
  		        		data[start+j]=D->shareF[n][m][iend-i][j]; 
  		      	start+=ny;
  		    	}  
  		MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
  	}             
	MPI_Barrier(MPI_COMM_WORLD);

 
	//Transferring odd ~ even cores 									 
	if(rank%2==0 && rank!=0)
	{
		MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
		start=0;
		for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
		  		for(i=1; i<share; i++) {
		    		for(j=0; j<ny; j++)
		      		D->shareF[n][m][istart-i][j]=data[start+j]; 
		      	start+=ny; 
		   	}
	}
	else if(rank%2==1 && rank!=D->L-1)
	{
	  	start=0; 
	  	for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)		  
	  	  		for(i=1; i<share; i++) {
	  	    		for(j=0; j<ny; j++)	
	  	      		data[start+j]=D->shareF[n][m][iend-i][j]; 
	  	      	start+=ny;
	  	   	}  
	  	MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
	}
	free(data);
}

void MPI_TransferFNew_Period_X(Domain *D,int numField,int ny,int share)
{
	int i,j,m,num,start,end,numMode,n;
	int istart,iend,jstart,jend;
	int myrank, nTasks,rank; 
	double *data;
	MPI_Status status;        	

	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);           	

	istart=D->istart;    iend=D->iend;
	jstart=D->jstart;    jend=D->jend;
	numMode=D->numMode;

	rank=myrank;	

	// Left to Right
	num=numMode*numField*ny*share;
	data = (double *)malloc(num*sizeof(double ));
 
	start=0; 
	for(n=0; n<numField; n++)
		for(m=0; m<numMode; m++)
			for(i=0; i<share; i++) {
				for(j=0; j<ny; j++) 
					data[start+j]=D->shareF[n][m][i+istart][j]; 
      		start+=ny; 
   		}
	if(D->L==1) {
		start=0; 
		for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
      		for(i=0; i<share; i++) {
        			for(j=0; j<ny; j++) 
						D->shareF[n][m][iend+i][j]=data[start+j]; 
          		start+=ny; 
      		}
	} else {
		if(rank==D->L-1)   {
      	MPI_Recv(data,num,MPI_DOUBLE,0,0, MPI_COMM_WORLD,&status);  
      	start=0; 
      	for(n=0; n<numField; n++)
				for(m=0; m<numMode; m++)
        			for(i=0; i<share; i++) {
          			for(j=0; j<ny; j++) 
              			D->shareF[n][m][iend+i][j]=data[start+j]; 
            		start+=ny; 
         		}
		}  else if(rank==0) {
      	MPI_Send(data,num,MPI_DOUBLE,D->L-1,myrank,MPI_COMM_WORLD);             
		}
   	MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);

  // Rgitht to Left 
	num=numMode*numField*ny*(share-1);
	data = (double *)malloc(num*sizeof(double ));

	start=0; 
	for(n=0; n<numField; n++)
		for(m=0; m<numMode; m++)	
   		for(i=1; i<share; i++) {
   			for(j=0; j<ny; j++)
         		data[start+j]=D->shareF[n][m][iend-i][j]; 
      		start+=ny; 
   		}

	if(D->L==1) {
		start=0;
   	for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)	
      		for(i=1; i<share; i++) {
      			for(j=0; j<ny; j++)
         			D->shareF[n][m][istart-i][j]=data[start+j]; 
         		start+=ny; 
        		}
	} else {      
		if(rank==0)   {
			MPI_Recv(data,num,MPI_DOUBLE,D->L-1,D->L-1, MPI_COMM_WORLD,&status);  
      	start=0;
      	for(n=0; n<numField; n++)
				for(m=0; m<numMode; m++)				
        			for(i=1; i<share; i++) {
          			for(j=0; j<ny; j++)
                  	D->shareF[n][m][istart-i][j]=data[start+j]; 
            		start+=ny; 
         		}
   	}  else if(rank==D->L-1)
      	MPI_Send(data,num,MPI_DOUBLE,0,myrank, MPI_COMM_WORLD);             
   	MPI_Barrier(MPI_COMM_WORLD);
	}
  
	free(data);
}



// ----------------------------------- Current, Density ----------------------------------------//

void MPI_TransferJNew_Xplus(Domain *D,int numField,int ny,int share)
{
   int i,j,m,start,numMode,num,n;
   int istart,iend,jstart,jend;
   int myrank, nTasks, rank;
   double *data;

   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;    iend=D->iend;
   jstart=D->jstart;    jend=D->jend;
   numMode=D->numMode;

   rank=myrank;
   num=numMode*numField*ny*share;
   data = (double *)malloc(num*sizeof(double ));

	//Transferring even ~ odd cores 
	if(rank%2==1)
	{
   	MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
   	start=0;
   	for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
      		for(i=0; i<share; i++) {
        			for(j=0; j<ny; j++) 
		        		D->shareF[n][m][istart+i][j]+=data[j+start]; 
			    	start+=ny; 
		   	} 
	}
	else if(rank%2==0 && rank!=D->L-1) 
	{
		start=0;
   	for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
      		for(i=0; i<share; i++) {
        			for(j=0; j<ny; j++)  {
         			data[start+j]=D->shareF[n][m][iend+i][j]; 
            		D->shareF[n][m][iend+i][j]=0.0;
         		} 
		   		start+=ny; 
		   	}
   	MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	//Transferring odd ~ even cores   
	if(rank%2==0 && rank!=0)
	{
		MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
		start=0;
		for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
      		for(i=0; i<share; i++) {
        			for(j=0; j<ny; j++) 
		      		D->shareF[n][m][istart+i][j]+=data[j+start]; 
			   	start+=ny; 
				} 
	}
	else if(rank%2==1 && rank!=D->L-1) 
	{
		start=0;
		for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
      		for(i=0; i<share; i++) {
        			for(j=0; j<ny; j++)  {
	  	      		data[start+j]=D->shareF[n][m][iend+i][j]; 
            		D->shareF[n][m][iend+i][j]=0.0;
         		} 
	  	   		start+=ny; 
	  			}
		MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	free(data);
}
                                                

void MPI_TransferJNew_Xminus(Domain *D,int numField,int ny,int share)
{
   int i,j,m,start,numMode,num,n;
   int istart,iend,jstart,jend;
   int myrank, nTasks, rank;
   double *data;

   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;    iend=D->iend;
   jstart=D->jstart;    jend=D->jend;
   numMode=D->numMode;

   rank=myrank;
   num=numMode*numField*ny*(share-1);
   data = (double *)malloc(num*sizeof(double ));

	//Transferring even ~ odd cores 
	if(rank%2==0 && rank!=D->L-1)
	{
		MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
   	start=0;
   	for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
      		for(i=1; i<share; i++) {
        			for(j=0; j<ny; j++) 
		        		D->shareF[n][m][iend-i][j]+=data[j+start]; 
			    	start+=ny; 
		   	} 
	}
	else if(rank%2==1) 
	{
   	start=0;
   	for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
      		for(i=1; i<share; i++) {
        			for(j=0; j<ny; j++)    {
	  	      		data[start+j]=D->shareF[n][m][istart-i][j]; 
            		D->shareF[n][m][istart-i][j]=0.0;
          		} 
	  	    		start+=ny; 
	  	  		}    
   	MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

   //Transferring odd ~ even cores             
	if(rank%2==1 && rank!=D->L-1)
	{
   	MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
   	start=0;
   	for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
      		for(i=1; i<share; i++) {
        			for(j=0; j<ny; j++) 
		        		D->shareF[n][m][iend-i][j]+=data[j+start]; 
			    	start+=ny; 
		   	} 
	}
	else if(rank%2==0 && rank!=0) 
	{
   	start=0;
   	for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
      		for(i=1; i<share; i++) {
        			for(j=0; j<ny; j++)   {
	  	      		data[start+j]=D->shareF[n][m][istart-i][j]; 
            		D->shareF[n][m][istart-i][j]=0.0; 
          		}
	  	    		start+=ny; 
	  			}    
   	MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	free(data);
}
                                                
                                        
void MPI_TransferJNew_Period_X(Domain *D,int numField,int ny,int share)
{
   int i,j,m,start,numMode,num,n;
   int istart,iend,jstart,jend;
   int myrank, nTasks, rank;
   double *data;

   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;    iend=D->iend;
   jstart=D->jstart;    jend=D->jend;
   numMode=D->numMode;

   rank=myrank;

	// Right to Left 
   num=numMode*numField*ny*share;
   data = (double *)malloc(num*sizeof(double ));

	start=0;
	for(n=0; n<numField; n++)
		for(m=0; m<numMode; m++)	
			for(i=0; i<share; i++) {
   			for(j=0; j<ny; j++)    {
      			data[start+j]=D->shareF[n][m][iend+i][j]; 
         		D->shareF[n][m][iend+i][j]=0.0; 
      		}
      		start+=ny; 
   		}

	if(D->L==1) {
   	start=0;
   	for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)	
      		for(i=0; i<share; i++) {
        			for(j=0; j<ny; j++) 
            		D->shareF[n][m][istart+i][j]+=data[j+start]; 
         		start+=ny; 
      		}			
	} else {
		if(rank==0)   {
      	MPI_Recv(data,num,MPI_DOUBLE,D->L-1,D->L-1,MPI_COMM_WORLD,&status);
      	start=0;
      	for(n=0; n<numField; n++)
				for(m=0; m<numMode; m++)
        			for(i=0; i<share; i++) {
          			for(j=0; j<ny; j++) 
              			D->shareF[n][m][istart+i][j]+=data[j+start]; 
            		start+=ny; 
         		}
   	}  else if(rank==D->L-1) {
      	MPI_Send(data,num,MPI_DOUBLE,0,myrank, MPI_COMM_WORLD);
		}
   	MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);

	// Left to Right 
   num=numMode*numField*ny*(share-1);
   data = (double *)malloc(num*sizeof(double ));

  	start=0;
	for(n=0; n<numField; n++)
		for(m=0; m<numMode; m++)
			for(i=1; i<share; i++) {
      		for(j=0; j<ny; j++) 	{
         		data[start+j]=D->shareF[n][m][istart-i][j]; 
         		D->shareF[n][m][istart-i][j]=0.0;
      		} 
      		start+=ny; 
   		}

	if(D->L==1) {
   	start=0;
   	for(n=0; n<numField; n++)
			for(m=0; m<numMode; m++)
      		for(i=1; i<share; i++) {
        			for(j=0; j<ny; j++) 
            		D->shareF[n][m][iend-i][j]+=data[j+start];
         		start+=ny; 
      		}
	} else {
		if(rank==D->L-1) {
      	MPI_Recv(data,num,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
      	start=0;
      	for(n=0; n<numField; n++)
				for(m=0; m<numMode; m++)
        			for(i=1; i<share; i++) {
          			for(j=0; j<ny; j++) 
              			D->shareF[n][m][iend-i][j]+=data[j+start]; 
         			start+=ny; 
         		}
   	} else if(rank==0) {
      	MPI_Send(data,num,MPI_DOUBLE,D->L-1,myrank, MPI_COMM_WORLD);
		}
   	MPI_Barrier(MPI_COMM_WORLD);
	}
   
	free(data);
}



void MPI_Transfer6F_Period_X(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    // Left to Right
    num=numMode*6*ny*share;
    data = (double *)malloc(num*sizeof(double ));
 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
      }
	
	if(nTasks==1)  {
		start=0; 
		for(m=0; m<numMode; m++)	  
			for(i=0; i<share; i++) {
				for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
		  		for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
			}
	} else {
		if(myrank==nTasks-1) {
			MPI_Recv(data,num,MPI_DOUBLE,0,0, MPI_COMM_WORLD,&status);  
		 	start=0; 
			for(m=0; m<numMode; m++)	  
				for(i=0; i<share; i++) {
					for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
		  			for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
				}
		} else if(myrank==0) {  
			MPI_Send(data,num,MPI_DOUBLE,nTasks-1,myrank,MPI_COMM_WORLD);
	 	}		  
    	MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);

    // Right to Left
    num=numMode*6*ny*(share-1);
    data = (double *)malloc(num*sizeof(double ));

	 start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
      }

	if(nTasks==1) {
		start=0;
		for(m=0; m<numMode; m++)
			for(i=1; i<share; i++) {
				for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
			}
	} else {
		if(myrank==0)   {
			MPI_Recv(data,num,MPI_DOUBLE,nTasks-1,nTasks-1, MPI_COMM_WORLD,&status);  
		  	start=0;
		  	for(m=0; m<numMode; m++)
				for(i=1; i<share; i++) {
					for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
				}
    	}
    	else if(myrank==nTasks-1)
      	MPI_Send(data,num,MPI_DOUBLE,0,myrank, MPI_COMM_WORLD);             
    	MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);
}

void MPI_Transfer8F_Period_X(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    // Left to Right
    num=numMode*8*ny*share;
    data = (double *)malloc(num*sizeof(double ));
	 //f1 ,f2, f3, f4, f5, f6, f7, f8, f9, f10,f11,f12
	 //ExR,ExI,BxR,BxI,PrR,PrI,PlR,PlI,SrR,SrI,SlR,SlI	 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][i+istart][j]; start+=ny;
      }
	
	if(nTasks==1)  {
		start=0; 
		for(m=0; m<numMode; m++)	  
			for(i=0; i<share; i++) {
				for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
		  		for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
			}
	} else {
		if(myrank==nTasks-1) {
			MPI_Recv(data,num,MPI_DOUBLE,0,0, MPI_COMM_WORLD,&status);  
		 	start=0; 
			for(m=0; m<numMode; m++)	  
				for(i=0; i<share; i++) {
					for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
		  			for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
				}
		} else if(myrank==0) {  
			MPI_Send(data,num,MPI_DOUBLE,nTasks-1,myrank,MPI_COMM_WORLD);
	 	}		  
    	MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);

    // Right to Left
    num=numMode*8*ny*(share-1);
    data = (double *)malloc(num*sizeof(double ));

	 //f1 ,f2, f3, f4, f5, f6, f7, f8, f9, f10,f11,f12
	 //ExR,ExI,BxR,BxI,PrR,PrI,PlR,PlI,SrR,SrI,SlR,SlI	 
	 start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
       for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][iend-i][j]; start+=ny;
      }

	if(nTasks==1) {
		start=0;
		for(m=0; m<numMode; m++)
			for(i=1; i<share; i++) {
				for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
			}
	} else {
		if(myrank==0)   {
			MPI_Recv(data,num,MPI_DOUBLE,nTasks-1,nTasks-1, MPI_COMM_WORLD,&status);  
		  	start=0;
		  	for(m=0; m<numMode; m++)
				for(i=1; i<share; i++) {
					for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
				}
    	}
    	else if(myrank==nTasks-1)
      	MPI_Send(data,num,MPI_DOUBLE,0,myrank, MPI_COMM_WORLD);             
    	MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);
}

/*
void MPI_Transfer8F_Period_X(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    // Left to Right
    num=numMode*8*ny*share;
    data = (double *)malloc(num*sizeof(double ));
	 //f1 ,f2, f3, f4, f5, f6, f7, f8, f9, f10,f11,f12
	 //ExR,ExI,BxR,BxI,PrR,PrI,PlR,PlI,SrR,SrI,SlR,SlI	 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=jstart+1; j<jend+3; j++) data[start+j-jstart-1]=f1[m][i+istart][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j-jstart-1]=f2[m][i+istart][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j-jstart-1]=f3[m][i+istart][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j-jstart-1]=f4[m][i+istart][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j-jstart-1]=f5[m][i+istart][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j-jstart-1]=f6[m][i+istart][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j-jstart-1]=f7[m][i+istart][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j-jstart-1]=f8[m][i+istart][j]; start+=ny;
      }
	
	if(nTasks==1)  {
		start=0; 
		for(m=0; m<numMode; m++)	  
			for(i=0; i<share; i++) {
				for(j=jstart+1; j<jend+3; j++) f1[m][iend+i][j]=data[start+j-jstart-1]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f2[m][iend+i][j]=data[start+j-jstart-1]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f3[m][iend+i][j]=data[start+j-jstart-1]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f4[m][iend+i][j]=data[start+j-jstart-1]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f5[m][iend+i][j]=data[start+j-jstart-1]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f6[m][iend+i][j]=data[start+j-jstart-1]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f7[m][iend+i][j]=data[start+j-jstart-1]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f8[m][iend+i][j]=data[start+j-jstart-1]; start+=ny;
			}
	} else {
		if(myrank==nTasks-1) {
			MPI_Recv(data,num,MPI_DOUBLE,0,0, MPI_COMM_WORLD,&status);  
		 	start=0; 
			for(m=0; m<numMode; m++)	  
				for(i=0; i<share; i++) {
					for(j=jstart+1; j<jend+3; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
				}
		} else if(myrank==0) {  
			MPI_Send(data,num,MPI_DOUBLE,nTasks-1,myrank,MPI_COMM_WORLD);
	 	}		  
    	MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);

    // Right to Left
    num=numMode*8*ny*(share-1);
    data = (double *)malloc(num*sizeof(double ));

	 //f1 ,f2, f3, f4, f5, f6, f7, f8, f9, f10,f11,f12
	 //ExR,ExI,BxR,BxI,PrR,PrI,PlR,PlI,SrR,SrI,SlR,SlI	 
	 start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=jstart+1; j<jend+3; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j]=f7[m][iend-i][j]; start+=ny;
        for(j=jstart+1; j<jend+3; j++) data[start+j]=f8[m][iend-i][j]; start+=ny;
      }

	if(nTasks==1) {
		start=0;
		for(m=0; m<numMode; m++)
			for(i=1; i<share; i++) {
				for(j=jstart+1; j<jend+3; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=jstart+1; j<jend+3; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
			}
	} else {
		if(myrank==0)   {
			MPI_Recv(data,num,MPI_DOUBLE,nTasks-1,nTasks-1, MPI_COMM_WORLD,&status);  
		  	start=0;
		  	for(m=0; m<numMode; m++)
				for(i=1; i<share; i++) {
					for(j=jstart+1; j<jend+3; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=jstart+1; j<jend+3; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
				}
    	}
    	else if(myrank==nTasks-1)
      	MPI_Send(data,num,MPI_DOUBLE,0,myrank, MPI_COMM_WORLD);             
    	MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);
}
*/
void MPI_Transfer12F_Period_X(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,double ***f9,double ***f10,double ***f11,double ***f12,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    // Left to Right
    num=numMode*12*ny*share;
    data = (double *)malloc(num*sizeof(double ));
	 //f1 ,f2, f3, f4, f5, f6, f7, f8, f9, f10,f11,f12
	 //ExR,ExI,BxR,BxI,PrR,PrI,PlR,PlI,SrR,SrI,SlR,SlI	 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f9[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f10[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f11[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f12[m][i+istart][j]; start+=ny;
      }
	
	if(nTasks==1)  {
		start=0; 
		for(m=0; m<numMode; m++)	  
			for(i=0; i<share; i++) {
				for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
		  		for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f9[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f10[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f11[m][iend+i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f12[m][iend+i][j]=data[start+j]; start+=ny;
			}
	} else {
		if(myrank==nTasks-1) {
			MPI_Recv(data,num,MPI_DOUBLE,0,0, MPI_COMM_WORLD,&status);  
		 	start=0; 
			for(m=0; m<numMode; m++)	  
				for(i=0; i<share; i++) {
					for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
		  			for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f9[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f10[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f11[m][iend+i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f12[m][iend+i][j]=data[start+j]; start+=ny;
				}
		} else if(myrank==0) {  
			MPI_Send(data,num,MPI_DOUBLE,nTasks-1,myrank,MPI_COMM_WORLD);
	 	}		  
    	MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);

    // Right to Left
    num=numMode*12*ny*(share-1);
    data = (double *)malloc(num*sizeof(double ));

	 //f1 ,f2, f3, f4, f5, f6, f7, f8, f9, f10,f11,f12
	 //ExR,ExI,BxR,BxI,PrR,PrI,PlR,PlI,SrR,SrI,SlR,SlI	 
	 start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
       for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f9[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f10[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f11[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f12[m][iend-i][j]; start+=ny;
      }

	if(nTasks==1) {
		start=0;
		for(m=0; m<numMode; m++)
			for(i=1; i<share; i++) {
				for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f9[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f10[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f11[m][istart-i][j]=data[start+j]; start+=ny;
				for(j=0; j<ny; j++) f12[m][istart-i][j]=data[start+j]; start+=ny;
			}
	} else {
		if(myrank==0)   {
			MPI_Recv(data,num,MPI_DOUBLE,nTasks-1,nTasks-1, MPI_COMM_WORLD,&status);  
		  	start=0;
		  	for(m=0; m<numMode; m++)
				for(i=1; i<share; i++) {
					for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f9[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f10[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f11[m][istart-i][j]=data[start+j]; start+=ny;
					for(j=0; j<ny; j++) f12[m][istart-i][j]=data[start+j]; start+=ny;
				}
    	}
    	else if(myrank==nTasks-1)
      	MPI_Send(data,num,MPI_DOUBLE,0,myrank, MPI_COMM_WORLD);             
    	MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);
}

void MPI_Transfer1F_Xminus(Domain *D,double ***f1,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*1*ny*3;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer1F_Xplus(Domain *D,double ***f1,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*1*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer2F_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*2*ny*3;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer2F_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*2*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer4F_Split_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*4*ny*(share+1);
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=-1; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=-1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=-1; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=-1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer4F_Split_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*4*ny*(share-2);
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=2; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=2; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=2; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=2; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}


void MPI_Transfer4F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*4*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer6F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*6*ny*3;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer6F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*6*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer8F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*8*ny*3;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 


    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==1) {
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f7[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f8[m][i+istart][j]; start+=ny;
        }
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);  
    }           
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             

    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==0 && rank!=0) {
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f7[m][i+istart][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f8[m][i+istart][j]; start+=ny;
        }
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD); 
    }            
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer8F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*8*ny*(share-1);
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 

      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1) {
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++)  {
          for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f7[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f8[m][iend-i][j]; start+=ny;
        }
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD); 
    }            
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores                     
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1) {
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++)  {
          for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f7[m][iend-i][j]; start+=ny;
          for(j=0; j<ny; j++) data[start+j]=f8[m][iend-i][j]; start+=ny;
        }
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);   
    }          
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer12F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,double ***f9,double ***f10,double ***f11,double ***f12,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*12*ny*share;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f9[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f10[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f11[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f12[m][i+istart][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f9[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f10[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f11[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f12[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f9[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f10[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f11[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f12[m][i+istart][j]; start+=ny;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f9[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f10[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f11[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f12[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer12F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,double ***f9,double ***f10,double ***f11,double ***f12,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*12*ny*(share-1);
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f9[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f10[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f11[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f12[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f9[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f10[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f11[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f12[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f9[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f10[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f11[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f12[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f9[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f10[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f11[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f12[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}


//--------------------------- Current -----------------------------//
/*
void MPI_TransferJ_Xplus(Domain *D,int numField)
{
  int i,j,m,n,start,numMode,num,ny,share;
  int istart,iend,jstart,jend;
  int myrank, nTasks, rank;
  double *data;

  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;
  ny=jend+3;
  share=3;

  rank=myrank/D->M;
  num=numField*numMode*ny*3;
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores 
  if(rank%2==1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
    start=0;
    for(n=0; n<numField; n++)
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) 
        {
          for(j=0; j<ny; j++) 
            D->shareF[n][m][istart+i][j]+=data[j+start];
          start+=ny;
        }
  }
  else if(rank%2==0 && rank!=D->L-1) {
    start=0;
    for(n=0; n<numField; n++)
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) 
        {
          for(j=0; j<ny; j++)  { 
            data[start+j]=D->shareF[n][m][iend+i][j];
            D->shareF[n][m][iend+i][j]=0.0;
          }
        }
    MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);

    
  //Transferring odd ~ even cores             
  if(rank%2==0 && rank!=0)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
    start=0;
    for(n=0; n<numField; n++)
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) 
        {
          for(j=0; j<ny; j++) 
            D->shareF[n][m][istart+i][j]+=data[j+start]; 
          start+=ny;
        }
  }
  else if(rank%2==1 && rank!=D->L-1) 
  {
    start=0;
    for(n=0; n<numField; n++)
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) 
        {
          for(j=0; j<ny; j++) {
            data[start+j]=D->shareF[n][m][istart+i][j];
            D->shareF[n][m][istart+i][j]=0.0;
          }
          start+=ny;
        }
    MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  free(data);
}


void MPI_TransferJ_Xminus(Domain *D,int numField)
{
  int i,j,n,m,start,numMode,num,ny,share;
  int istart,iend,jstart,jend;
  int myrank, nTasks, rank;
  double *data;

  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;
  ny=jend+3;
  share=3;

  rank=myrank/D->M;
  num=numField*numMode*ny*(share-1);
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores 
  if(rank%2==0 && rank!=D->L-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
    start=0;
    for(n=0; n<numField; n++)
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) 
        {
          for(j=0; j<ny; j++) 
            D->shareF[n][m][iend-i][j]+=data[j+start]; 
          start+=ny;
        }
  }
  else if(rank%2==1)   {
    start=0;
    for(n=0; n<numField; n++)
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) 
        {
          for(j=0; j<ny; j++) { 
            data[start+j]=D->shareF[n][m][istart-i][j]; 
            D->shareF[n][m][istart-i][j]=0.0; 
          } 
          start+=ny;
        }
    MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);


  //Transferring odd ~ even cores             
  if(rank%2==1 && rank!=D->L-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
    start=0;
    for(n=0; n<numField; n++)
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) 
        {
          for(j=0; j<ny; j++) 
            D->shareF[n][m][iend-i][j]+=data[j+start]; 
          start+=ny;
        }
  }
  else if(rank%2==0 && rank!=0) 
  {
    start=0;
    for(n=0; n<numField; n++)
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) 
        {
          for(j=0; j<ny; j++) { 
            data[start+j]=D->shareF[n][m][istart-i][j]; 
            D->shareF[n][m][istart-i][j]=0.0; 
          } 
          start+=ny;
        }
    MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  free(data);
}

void MPI_TransferJ_Period_X(Domain *D,int numField)
{
  int i,j,m,n,start,numMode,num,ny,share;
  int istart,iend,jstart,jend;
  int myrank, nTasks, rank;
  double *data;

  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;
  ny=jend+3;
  share=3;

  rank=myrank/D->M;
  num=numField*numMode*ny*share;
  data = (double *)malloc(num*sizeof(double ));

  // Right to Left
  start=0;
  for(n=0; n<numField; n++)
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) 
      {
        for(j=0; j<ny; j++) { 
          data[start+j]=D->shareF[n][m][iend+i][j]; 
          D->shareF[n][m][iend+i][j]=0.0; 
        } 
        start+=ny;
      }

	if(nTasks==1) {
   	start=0;
    for(n=0; n<numField; n++)
      for(m=0; m<numMode; m++)
      	for(i=0; i<share; i++) 
        {
         	for(j=0; j<ny; j++) 
            D->shareF[n][m][istart+i][j]+=data[j+start]; 
          start+=ny;
        }
	} else {
		if(myrank==0)
		{
      MPI_Recv(data,num,MPI_DOUBLE,nTasks-1,nTasks-1,MPI_COMM_WORLD,&status);
      start=0;
      for(n=0; n<numField; n++)
        for(m=0; m<numMode; m++)
       	  for(i=0; i<share; i++) 
          {
         		for(j=0; j<ny; j++) 
              D->shareF[n][m][istart+i][j]+=data[j+start]; 
            start+=ny;
       	  }
    }
		else if(myrank==nTasks-1)
      MPI_Send(data,num,MPI_DOUBLE,0,myrank, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);

  // Left to Right
  num=numField*numMode*ny*(share-1);
  data = (double *)malloc(num*sizeof(double ));

  start=0;
  for(n=0; n<numField; n++)
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++) 
      {
        for(j=0; j<ny; j++) { 
          data[start+j]=D->shareF[n][m][istart-i][j]; 
          D->shareF[n][m][istart-i][j]=0.0; 
        } 
        start+=ny;
      }

	if(nTasks==1) {
    start=0;
    for(n=0; n<numField; n++)
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) 
        {
          for(j=0; j<ny; j++) 
            D->shareF[n][m][iend-i][j]+=data[j+start]; 
          start+=ny;
        }
	} else {
		if(myrank==nTasks-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
      start=0;
      for(n=0; n<numField; n++)
			  for(m=0; m<numMode; m++)
         	for(i=1; i<share; i++) 
          {
           	for(j=0; j<ny; j++) 
              D->shareF[n][m][iend-i][j]+=data[j+start]; 
            start+=ny;
      		}
    }
		else if(myrank==0)
      MPI_Send(data,num,MPI_DOUBLE,nTasks-1,myrank, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
	}

	free(data);
}
*/


void MPI_TransferJ_Period_X(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share)
{
    int i,j,m,start,numMode,num;
    int istart,iend,jstart,jend;
    int myrank, nTasks, rank;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    rank=myrank/D->M;
    num=numMode*6*ny*3;
    data = (double *)malloc(num*sizeof(double ));

    // Right to Left
    start=0;
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) { data[start+j]=f1[m][iend+i][j]; f1[m][iend+i][j]=0.0; } start+=ny;
        for(j=0; j<ny; j++) { data[start+j]=f2[m][iend+i][j]; f2[m][iend+i][j]=0.0; } start+=ny;
        for(j=0; j<ny; j++) { data[start+j]=f3[m][iend+i][j]; f3[m][iend+i][j]=0.0; } start+=ny;
        for(j=0; j<ny; j++) { data[start+j]=f4[m][iend+i][j]; f4[m][iend+i][j]=0.0; } start+=ny;
        for(j=0; j<ny; j++) { data[start+j]=f5[m][iend+i][j]; f5[m][iend+i][j]=0.0; } start+=ny;
        for(j=0; j<ny; j++) { data[start+j]=f6[m][iend+i][j]; f6[m][iend+i][j]=0.0; } start+=ny;
      }

	if(nTasks==1) {
   	start=0;
      for(m=0; m<numMode; m++)
      	for(i=0; i<share; i++) {
         	for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           	for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
           	for(j=0; j<ny; j++) f3[m][istart+i][j]+=data[j+start]; start+=ny;
           	for(j=0; j<ny; j++) f4[m][istart+i][j]+=data[j+start]; start+=ny;
           	for(j=0; j<ny; j++) f5[m][istart+i][j]+=data[j+start]; start+=ny;
           	for(j=0; j<ny; j++) f6[m][istart+i][j]+=data[j+start]; start+=ny;
        	}
	} else {
		if(myrank==0)
		{
      	MPI_Recv(data,num,MPI_DOUBLE,nTasks-1,nTasks-1,MPI_COMM_WORLD,&status);
       	start=0;
       	for(m=0; m<numMode; m++)
         	for(i=0; i<share; i++) {
           		for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           		for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
		         for(j=0; j<ny; j++) f3[m][istart+i][j]+=data[j+start]; start+=ny;
      		   for(j=0; j<ny; j++) f4[m][istart+i][j]+=data[j+start]; start+=ny;
		         for(j=0; j<ny; j++) f5[m][istart+i][j]+=data[j+start]; start+=ny;
      		   for(j=0; j<ny; j++) f6[m][istart+i][j]+=data[j+start]; start+=ny;
         	}
    	}
		else if(myrank==nTasks-1)
      	MPI_Send(data,num,MPI_DOUBLE,0,myrank, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);

    // Left to Right
    num=numMode*6*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    start=0;
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++) {
        for(j=0; j<ny; j++) { data[start+j]=f1[m][istart-i][j]; f1[m][istart-i][j]=0.0; } start+=ny;
        for(j=0; j<ny; j++) { data[start+j]=f2[m][istart-i][j]; f2[m][istart-i][j]=0.0; } start+=ny;
        for(j=0; j<ny; j++) { data[start+j]=f3[m][istart-i][j]; f3[m][istart-i][j]=0.0; } start+=ny;
        for(j=0; j<ny; j++) { data[start+j]=f4[m][istart-i][j]; f4[m][istart-i][j]=0.0; } start+=ny;
        for(j=0; j<ny; j++) { data[start+j]=f5[m][istart-i][j]; f5[m][istart-i][j]=0.0; } start+=ny;
        for(j=0; j<ny; j++) { data[start+j]=f6[m][istart-i][j]; f6[m][istart-i][j]=0.0; } start+=ny;
      }

	if(nTasks==1) {
       start=0;
       for(m=0; m<numMode; m++)
         for(i=1; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f3[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f4[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f5[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f6[m][iend-i][j]+=data[j+start]; start+=ny;
         }
	} else {
		if(myrank==nTasks-1)
    	{
      	MPI_Recv(data,num,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
       	start=0;
			for(m=0; m<numMode; m++)
         	for(i=1; i<share; i++) {
           		for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
		         for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
      		   for(j=0; j<ny; j++) f3[m][iend-i][j]+=data[j+start]; start+=ny;
		         for(j=0; j<ny; j++) f4[m][iend-i][j]+=data[j+start]; start+=ny;
      		   for(j=0; j<ny; j++) f5[m][iend-i][j]+=data[j+start]; start+=ny;
		         for(j=0; j<ny; j++) f6[m][iend-i][j]+=data[j+start]; start+=ny;
      		}
    	}
		else if(myrank==0)
      	MPI_Send(data,num,MPI_DOUBLE,nTasks-1,myrank, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
	}

	free(data);
}

void MPI_TransferJ_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share)
{
    int i,j,m,start,numMode,num;
    int istart,iend,jstart,jend;
    int myrank, nTasks, rank;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    rank=myrank/D->M;
    num=numMode*6*ny*3;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 


    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=0; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f3[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f4[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f5[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f6[m][istart+i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==0 && rank!=D->L-1) {
      start=0;
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) { data[start+j]=f1[m][iend+i][j]; f1[m][iend+i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f2[m][iend+i][j]; f2[m][iend+i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f3[m][iend+i][j]; f3[m][iend+i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f4[m][iend+i][j]; f4[m][iend+i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f5[m][iend+i][j]; f5[m][iend+i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f6[m][iend+i][j]; f6[m][iend+i][j]=0.0; } start+=ny;
        }
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             


    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=0; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f3[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f4[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f5[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f6[m][istart+i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==1 && rank!=D->L-1) {
      start=0;
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) { data[start+j]=f1[m][iend+i][j]; f1[m][iend+i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f2[m][iend+i][j]; f2[m][iend+i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f3[m][iend+i][j]; f3[m][iend+i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f4[m][iend+i][j]; f4[m][iend+i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f5[m][iend+i][j]; f5[m][iend+i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f6[m][iend+i][j]; f6[m][iend+i][j]=0.0; } start+=ny;
        }
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    free(data);
}
                                                
void MPI_TransferJ_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share)
{
    int i,j,m,start,numMode,num;
    int istart,iend,jstart,jend;
    int myrank, nTasks, rank;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    rank=myrank/D->M;
    num=numMode*6*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    if(rank%2==0 && rank!=D->L-1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=1; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f3[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f4[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f5[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f6[m][iend-i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==1)   {
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) { data[start+j]=f1[m][istart-i][j]; f1[m][istart-i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f2[m][istart-i][j]; f2[m][istart-i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f3[m][istart-i][j]; f3[m][istart-i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f4[m][istart-i][j]; f4[m][istart-i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f5[m][istart-i][j]; f5[m][istart-i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f6[m][istart-i][j]; f6[m][istart-i][j]=0.0; } start+=ny;
        }
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);


    //Transferring odd ~ even cores             
    if(rank%2==1 && rank!=D->L-1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=1; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f3[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f4[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f5[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f6[m][iend-i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==0 && rank!=0) {
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) { data[start+j]=f1[m][istart-i][j]; f1[m][istart-i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f2[m][istart-i][j]; f2[m][istart-i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f3[m][istart-i][j]; f3[m][istart-i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f4[m][istart-i][j]; f4[m][istart-i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f5[m][istart-i][j]; f5[m][istart-i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f6[m][istart-i][j]; f6[m][istart-i][j]=0.0; } start+=ny;
        }
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
   
    free(data);
}

// void MPI_TransferRho_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share)
// {
//     int i,j,m,start,numMode,num;
//     int istart,iend,jstart,jend;
//     int myrank, nTasks, rank;
//     double *data;

//     MPI_Status status;

//     MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
//     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

//     istart=D->istart;    iend=D->iend;
//     jstart=D->jstart;    jend=D->jend;
//     numMode=D->numMode;

//     rank=myrank/D->M;
//     num=numMode*2*ny*3;
//     data = (double *)malloc(num*sizeof(double ));

//     //Transferring even ~ odd cores 
//     start=0;
//     for(m=0; m<numMode; m++)
//       for(i=0; i<share; i++) {
//         for(j=0; j<ny; j++) data[start+j]=f1[m][iend+i][j]; start+=ny;
//         for(j=0; j<ny; j++) data[start+j]=f2[m][iend+i][j]; start+=ny;
//       }

//     if(rank%2==1)
//     {
//        MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
//        start=0;
//        for(m=0; m<numMode; m++)
//          for(i=0; i<share; i++) {
//            for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
//            for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
//          }
//     }
//     else if(rank%2==0 && rank!=D->L-1)
//       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);

//     MPI_Barrier(MPI_COMM_WORLD);

//     //Transferring odd ~ even cores             
//     start=0;
//     for(m=0; m<numMode; m++)
//       for(i=0; i<share; i++) {
//         for(j=0; j<ny; j++) data[start+j]=f1[m][iend+i][j]; start+=ny;
//         for(j=0; j<ny; j++) data[start+j]=f2[m][iend+i][j]; start+=ny;
//       }

//     if(rank%2==0 && rank!=0)
//     {
//        MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
//        start=0;
//        for(m=0; m<numMode; m++)
//          for(i=0; i<share; i++) {
//            for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
//            for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
//          }
//     }
//     else if(rank%2==1 && rank!=D->L-1)
//       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
//     MPI_Barrier(MPI_COMM_WORLD);

//     free(data);
// }
                                                
// void MPI_TransferRho_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share)
// {
//     int i,j,m,start,numMode,num;
//     int istart,iend,jstart,jend;
//     int myrank, nTasks, rank;
//     double *data;

//     MPI_Status status;

//     MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
//     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

//     istart=D->istart;    iend=D->iend;
//     jstart=D->jstart;    jend=D->jend;
//     numMode=D->numMode;

//     rank=myrank/D->M;
//     num=numMode*2*ny*2;
//     data = (double *)malloc(num*sizeof(double ));

//     //Transferring even ~ odd cores 
//     start=0;
//     for(m=0; m<numMode; m++)
//       for(i=1; i<share; i++) {
//         for(j=0; j<ny; j++) data[start+j]=f1[m][istart-i][j]; start+=ny;
//         for(j=0; j<ny; j++) data[start+j]=f2[m][istart-i][j]; start+=ny;
//       }

//     if(rank%2==0 && rank!=D->L-1)
//     {
//        MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
//        start=0;
//        for(m=0; m<numMode; m++)
//          for(i=1; i<share; i++) {
//            for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
//            for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
//          }
//     }
//     else if(rank%2==1)
//       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);

//     MPI_Barrier(MPI_COMM_WORLD);

//     //Transferring odd ~ even cores             
//     start=0;
//     for(m=0; m<numMode; m++)
//       for(i=1; i<share; i++) {
//         for(j=0; j<ny; j++) data[start+j]=f1[m][istart-i][j]; start+=ny;
//         for(j=0; j<ny; j++) data[start+j]=f2[m][istart-i][j]; start+=ny;
//       }

//     if(rank%2==1 && rank!=D->L-1)
//     {
//        MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
//        start=0;
//        for(m=0; m<numMode; m++)
//          for(i=1; i<share; i++) {
//            for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
//            for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
//          }
//     }
//     else if(rank%2==0 && rank!=0)
//       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
//     MPI_Barrier(MPI_COMM_WORLD);
   
//     free(data);
// }

void MPI_TransferDen_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share)
{
    int i,j,m,start,numMode,num;
    int istart,iend,jstart,jend;
    int myrank, nTasks, rank;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    rank=myrank/D->M;
    num=numMode*2*ny*3;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 


    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=0; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==0 && rank!=D->L-1) {
      start=0;
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) { data[start+j]=f1[m][iend+i][j]; f1[m][iend+i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f2[m][iend+i][j]; f2[m][iend+i][j]=0.0; } start+=ny;
        }
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=0; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==1 && rank!=D->L-1) {
      start=0;
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) { data[start+j]=f1[m][iend+i][j]; f1[m][iend+i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f2[m][iend+i][j]; f2[m][iend+i][j]=0.0; } start+=ny;
        }
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    }
 
    MPI_Barrier(MPI_COMM_WORLD);
 
    free(data);
}
                                              
void MPI_TransferDen_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share)
{
    int i,j,m,start,numMode,num;
    int istart,iend,jstart,jend;
    int myrank, nTasks, rank;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    rank=myrank/D->M;
    num=numMode*2*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 


    if(rank%2==0 && rank!=D->L-1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=1; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==1) {
    start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) { data[start+j]=f1[m][istart-i][j]; f1[m][istart-i][j]=0.0; } start+=ny;
          for(j=0; j<ny; j++) { data[start+j]=f2[m][istart-i][j]; f2[m][istart-i][j]=0.0; } start+=ny;
        }
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
    }


    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    if(rank%2==1 && rank!=D->L-1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=1; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==0 && rank!=0) {
      start=0;
        for(m=0; m<numMode; m++)
          for(i=1; i<share; i++) {
            for(j=0; j<ny; j++) { data[start+j]=f1[m][istart-i][j]; f1[m][istart-i][j]=0.0; } start+=ny;
            for(j=0; j<ny; j++) { data[start+j]=f2[m][istart-i][j]; f2[m][istart-i][j]=0.0; } start+=ny;
          }
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
   
    free(data);
}

void MPI_TransferDen_Period_X(Domain *D,double ***f1,double ***f2,int ny,int share)                                                
{
    int i,j,m,start,numMode,num;
    int istart,iend,jstart,jend;
    int myrank, nTasks, rank;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    rank=myrank/D->M;
    num=numMode*2*ny*3;
    data = (double *)malloc(num*sizeof(double ));

    // Right to Left
    start=0;
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) { data[start+j]=f1[m][iend+i][j]; f1[m][iend+i][j]=0.0; } start+=ny;
        for(j=0; j<ny; j++) { data[start+j]=f2[m][iend+i][j]; f2[m][iend+i][j]=0.0; } start+=ny;
      }

	if(nTasks==1) {
   	start=0;
      for(m=0; m<numMode; m++)
      	for(i=0; i<share; i++) {
         	for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
        }
	} else {
		if(myrank==0)
		{
      	MPI_Recv(data,num,MPI_DOUBLE,nTasks-1,nTasks-1,MPI_COMM_WORLD,&status);
       	start=0;
       	for(m=0; m<numMode; m++)
         	for(i=0; i<share; i++) {
           		for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           		for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
         	}
    }
		else if(myrank==nTasks-1)
      MPI_Send(data,num,MPI_DOUBLE,0,myrank, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
  }
	free(data);

  // Left to Right
  num=numMode*2*ny*2;
  data = (double *)malloc(num*sizeof(double ));

  start=0;
  for(m=0; m<numMode; m++)
    for(i=1; i<share; i++) {
      for(j=0; j<ny; j++) { data[start+j]=f1[m][istart-i][j]; f1[m][istart-i][j]=0.0; } start+=ny;
      for(j=0; j<ny; j++) { data[start+j]=f2[m][istart-i][j]; f2[m][istart-i][j]=0.0; } start+=ny;
    }

	if(nTasks==1) {
       start=0;
       for(m=0; m<numMode; m++)
         for(i=1; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
         }
	} else {
		if(myrank==nTasks-1)
    	{
      	MPI_Recv(data,num,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
       	start=0;
			for(m=0; m<numMode; m++)
         	for(i=1; i<share; i++) {
           		for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
		         for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
      		}
    	}
		else if(myrank==0)
      	MPI_Send(data,num,MPI_DOUBLE,nTasks-1,myrank, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
	}

	free(data);
}