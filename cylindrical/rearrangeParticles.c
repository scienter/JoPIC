#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"


void rearrangeParticles(Domain *D)
{
    Particle **particle;
    particle=D->particle;

    int i,j,s,intX=0,intY=0,cnt,deleteFlag=0;
    int istart,iend,jstart,jend;
    double z,x,y,R,r;
    ptclList *p,*New,*prev,*tmp;
    int myrank;
  	 MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<D->nSpecies; s++)
        {
          cnt=1;
          p=particle[i][j].head[s]->pt;
          while(p)
          {
            if(cnt==1)
              prev=p;
            else	;
            deleteFlag=0;
              
            z=p->z; x=p->x; y=p->y; R=sqrt(x*x+y*y);
            r=R-(j-jstart);

            if(z>=1.0)  {
              intX=(int)z;
              z-=intX;
              deleteFlag=1;
            } else if(z<0.0) {              
              intX=(int)(z-1.0);
              z-=intX;
              deleteFlag=1;
            } else  intX=0;

            if(r>=1.0)  {
              intY=(int)r;
              deleteFlag=1;
            } else if(r<0.0) {              
              intY=(int)(r-1.0);
              deleteFlag=1;
            } else  intY=0;

            if(deleteFlag==1)
            {
              if(cnt==1)
              {
                p->z=z;
                particle[i][j].head[s]->pt = p->next;
                p->next = particle[i+intX][j+intY].head[s]->pt;
                particle[i+intX][j+intY].head[s]->pt = p;
                p=particle[i][j].head[s]->pt;
                cnt=1;
              }
              else
              {
                p->z=z;
                prev->next = p->next;
                p->next = particle[i+intX][j+intY].head[s]->pt;
                particle[i+intX][j+intY].head[s]->pt = p;
                p=prev->next;
              }
            }		//End of if(deleteFlag==1)
            else
            {
              prev=p;
              p=p->next;
              cnt++;
            }              
          }		//End if while(p)
        }		//End of for(s)
}

void rearrangeParticles_Period_Y(Domain *D)
{
  Particle **particle;
  particle=D->particle;

  int i,j,s,intX=0,intY=0,cnt,deleteFlag=0;
  int istart,iend,jstart,jend;
  double z,x,y,R,r,sinTh,cosTh;
  ptclList *p,*New,*prev,*tmp;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;

  j=jend;
  for(i=istart; i<iend; i++)
    for(s=0; s<D->nSpecies; s++)
    {
      cnt=1;
      p=particle[i][j].head[s]->pt;
      while(p)
      {
        if(cnt==1)
          prev=p;
        else	;
           
        z=p->z; x=p->x; y=p->y; R=sqrt(x*x+y*y);
        cosTh=x/R; sinTh=y/R;
        r=R-(j-jstart);
        R-=r;
        p->x=cosTh*R; p->y=sinTh*R;
        p->px*=-1; p->py*=-1;

        if(cnt==1)
        {
          particle[i][j].head[s]->pt = p->next;
          p->next = particle[i][j-1].head[s]->pt;
          particle[i][j-1].head[s]->pt = p;
          p=particle[i][j].head[s]->pt;
          cnt=1;
        }
        else
        {
          prev->next = p->next;
          p->next = particle[i][j-1].head[s]->pt;
          particle[i][j-1].head[s]->pt = p;
          p=prev->next;
        }
      }
    }
}
