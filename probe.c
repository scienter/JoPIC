#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

void trackID(Domain *D,int iteration,int istart,int iend,int jstart,int jend,int kstart,int kend)
{
  int i,n,iter,s,ii,index;
  ptclList *p;
  int species[D->idNums],ids[D->idNums],cores[D->idNums];

  for(i=0; i<D->idNums; i++)
  {
    ids[i]=D->trackID[i];
    cores[i]=D->trackCore[i];
    species[i]=D->trackS[i];
  }

  for(n=0; n<D->idNums; n++)
  {
    s=species[n];
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {
            if(p->id==ids[n] && p->core==cores[n])
            {
              D->track[n][iteration].x=p->x;
              D->track[n][iteration].y=p->y;
              D->track[n][iteration].z=p->z;
              D->track[n][iteration].px=p->px;
              D->track[n][iteration].py=p->py;
              D->track[n][iteration].pz=p->pz;
              D->track[n][iteration].id=p->id;
              D->track[n][iteration].core=p->core;
            }
            p=p->next;
          }
        }
  }
}

