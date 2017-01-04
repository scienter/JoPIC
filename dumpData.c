#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void saveDump(Domain D,int iteration)
{
  void saveDumpData();

  switch ((D.fieldType-1)*3+D.dimension)  {
  case ((1-1)*3+2) :
    saveDumpData(&D,iteration,D.istart,D.iend,D.jstart,D.jend,0,1,D.nxSub+5,D.nySub+5,1);
    break;
  case ((1-1)*3+3) :
    saveDumpData(&D,iteration,D.istart,D.iend,D.jstart,D.jend,D.kstart,D.kend,D.nxSub+5,D.nySub+5,D.nzSub+5);
    break;
  }

}

void restoreData(Domain *D,int iteration)
{
   void restore();

  switch ((D->fieldType-1)*3+D->dimension)  {
  case ((1-1)*3+2) :
    restore(D,iteration,D->istart,D->iend,D->jstart,D->jend,0,1,D->nxSub+5,D->nySub+5,1);
    break;
  case ((1-1)*3+3) :
    restore(D,iteration,D->istart,D->iend,D->jstart,D->jend,D->kstart,D->kend,D->nxSub+5,D->nySub+5,D->nzSub+5);
    break;
  }

}

void saveDumpData(Domain *D,int iteration,int istart,int iend,int jstart,int jend,int kstart,int kend,int totalX,int totalY,int totalZ)
{
   FILE *out;
   char name[100];
   int i,j,k,n,s,cnt;
   Particle ***particle;
   particle=D->particle;
   LoadList *LL;
//   Probe **probe;
//   probe=D.probe;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   sprintf(name,"dump%d_%d",iteration,myrank);
   out = fopen(name, "w");   

   // Save simulation Domain information
   fwrite(&(D->minXSub),sizeof(int),1,out);
   fwrite(&(D->maxXSub),sizeof(int),1,out);
   fwrite(&(D->minYSub),sizeof(int),1,out);
   fwrite(&(D->maxYSub),sizeof(int),1,out);
   fwrite(&(D->minZSub),sizeof(int),1,out);
   fwrite(&(D->maxZSub),sizeof(int),1,out);
   fwrite(&(D->nSpecies),sizeof(int),1,out);   
   LL=D->loadList;
   for(i=0; i<D->nSpecies; i++)
   {
     fwrite(&(LL->index),sizeof(int),1,out);   
     LL=LL->next;
   }
     
   // Save informations of particles inside the domain
   for (i=istart; i<iend; i++) 
     for (j=jstart; j<jend; j++)
       for (k=kstart; k<kend; k++)
         for(s=0; s<D->nSpecies; s++)
         {
           p = particle[i][j][k].head[s]->pt;
           cnt = 0;
           while(p)  { 
             p=p->next; 
             cnt++;    
           }
           fwrite(&cnt,sizeof(int),1,out);

           p = particle[i][j][k].head[s]->pt;
           while(p)  
           { 
             fwrite(&(p->x),sizeof(double),1,out);   
             fwrite(&(p->y),sizeof(double),1,out);   
             fwrite(&(p->z),sizeof(double),1,out);   
             fwrite(&(p->p1),sizeof(double),1,out);   
             fwrite(&(p->p2),sizeof(double),1,out);   
             fwrite(&(p->p3),sizeof(double),1,out);   
             fwrite(&(p->index),sizeof(int),1,out);   
             fwrite(&(p->core),sizeof(int),1,out);   

             p = p->next; 
           }            
         }		//End of for(s,i,j)

     for(i=0; i<totalX; i++) 
       for(j=0; j<totalY; j++)
         for(k=0; k<totalZ; k++)
         {
           fwrite(&(D->Ex[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->Pr[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->Pl[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->Bx[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->Sr[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->Sl[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->ExC[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->PrC[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->PlC[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->BxC[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->SrC[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->SlC[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->Jx[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->Jy[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->Jz[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->JxOld[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->JyOld[i][j][k]),sizeof(double),1,out);   
           fwrite(&(D->JzOld[i][j][k]),sizeof(double),1,out);   
         }

/*  
   //save Probe data
   if(D.probeNum>0)
   {
     for(n=0; n<D.probeNum; n++)
       for(i=0; i<=D.maxStep; i++)
       { 
         fwrite(&(probe[n][i].Pr),sizeof(double),1,out);   
         fwrite(&(probe[n][i].Pl),sizeof(double),1,out);   
         fwrite(&(probe[n][i].Sr),sizeof(double),1,out);   
         fwrite(&(probe[n][i].Sl),sizeof(double),1,out);   
         fwrite(&(probe[n][i].E1),sizeof(double),1,out);   
         fwrite(&(probe[n][i].B1),sizeof(double),1,out);   
       }
   }
*/
   fclose(out);
}

// Temporary routine to dump and restore data
void restore(Domain *D, int iteration,int istart,int iend,int jstart,int jend,int kstart,int kend,int totalX,int totalY,int totalZ)
{
   FILE *in;
   char name[100];
   int i,j,k,n,s,cnt,maxStep,nSpecies;
   double tmp;
   ptclList *p;
   Particle ***particle;
   particle=D->particle;
   LoadList *LL;
   Probe **probe;
   probe=D->probe;
   int myrank, nTasks;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   sprintf(name,"dump%d_%d",iteration,myrank);
   in = fopen(name, "r");   

   // restore simulation Domain information
    fread(&(D->minXSub),sizeof(int),1,in);
    fread(&(D->maxXSub),sizeof(int),1,in);
    fread(&(D->minYSub),sizeof(int),1,in);
    fread(&(D->maxYSub),sizeof(int),1,in);
    fread(&(D->minZSub),sizeof(int),1,in);
    fread(&(D->maxZSub),sizeof(int),1,in);
    fread(&(nSpecies),sizeof(int),1,in);
    LL=D->loadList;
    for(i=0; i<nSpecies; i++)
    {
      fread(&(LL->index),sizeof(int),1,in);
      LL=LL->next;
    }
//    D->minXSub+=1;
//    D->maxXSub+=1;

   // restore informations of particles inside the domain
    for (i=istart; i<iend; i++) 
      for (j=jstart; j<jend; j++)
        for (k=kstart; k<kend; k++)
        for(s=0; s<nSpecies; s++)
        {
          fread(&cnt,sizeof(int),1,in);

          for(n=0; n<cnt; n++)  
          { 
            p = (ptclList *)malloc(sizeof(ptclList)); 
            p->next = particle[i][j][k].head[s]->pt;
            particle[i][j][k].head[s]->pt = p;

            fread(&(p->x),sizeof(double),1,in);   
            fread(&(p->y),sizeof(double),1,in);   
            fread(&(p->z),sizeof(double),1,in);   
            fread(&(p->p1),sizeof(double),1,in);   
            fread(&(p->p2),sizeof(double),1,in);   
            fread(&(p->p3),sizeof(double),1,in);   
            fread(&(p->index),sizeof(int),1,in); 
            fread(&(p->core),sizeof(int),1,in); 
          }
        }

    for(i=0; i<totalX; i++) 
      for(j=0; j<totalY; j++)
        for(k=0; k<totalZ; k++)
        {
          fread(&(D->Ex[i][j][k]),sizeof(double),1,in);
          fread(&(D->Pr[i][j][k]),sizeof(double),1,in);
          fread(&(D->Pl[i][j][k]),sizeof(double),1,in);
          fread(&(D->Bx[i][j][k]),sizeof(double),1,in);
          fread(&(D->Sr[i][j][k]),sizeof(double),1,in);
          fread(&(D->Sl[i][j][k]),sizeof(double),1,in);
          fread(&(D->ExC[i][j][k]),sizeof(double),1,in);
          fread(&(D->PrC[i][j][k]),sizeof(double),1,in);
          fread(&(D->PlC[i][j][k]),sizeof(double),1,in);
          fread(&(D->BxC[i][j][k]),sizeof(double),1,in);
          fread(&(D->SrC[i][j][k]),sizeof(double),1,in);
          fread(&(D->SlC[i][j][k]),sizeof(double),1,in);
          fread(&(D->Jx[i][j][k]),sizeof(double),1,in);
          fread(&(D->Jy[i][j][k]),sizeof(double),1,in);
          fread(&(D->Jz[i][j][k]),sizeof(double),1,in);
          fread(&(D->JxOld[i][j][k]),sizeof(double),1,in);
          fread(&(D->JyOld[i][j][k]),sizeof(double),1,in);
          fread(&(D->JzOld[i][j][k]),sizeof(double),1,in);
        }

/*
   //restore Probe data
   if(D->probeNum>0)
   {
     for(n=0; n<D->probeNum; n++)
       for(i=0; i<=maxStep; i++)
       { 
         fread(&(probe[n][i].Pr),sizeof(double),1,in);   
         fread(&(probe[n][i].Pl),sizeof(double),1,in);   
         fread(&(probe[n][i].Sr),sizeof(double),1,in);   
         fread(&(probe[n][i].Sl),sizeof(double),1,in);   
         fread(&(probe[n][i].E1),sizeof(double),1,in);   
         fread(&(probe[n][i].B1),sizeof(double),1,in);   
       }
   }
*/   
    fclose(in);
}

