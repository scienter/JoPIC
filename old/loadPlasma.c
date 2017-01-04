#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>

float maxwellianVelocity(float temperature);
float randomValue(float beta);
void loadPlasma_crystal(Domain *D,LoadList *LL,int s);

float channel(float centerY,float y,float centerZ,float z,float coef)
{
  float result;
  result=1.0+coef*((y-centerY)*(y-centerY)+(z-centerZ)*(z-centerZ));
  return result;
}

void loadPlasma(Domain *D,LoadList *LL,int s,int iteration)
{
  void loadPolygonPlasma2D();
  void loadPolygonPlasma3D();
  void loadChannelPlasma();

  switch((LL->type-1)*3+D->dimension)
  {
  case ((Polygon-1)*3+1):
    loadPlasma_crystal(D,LL,s);
    break;
  case ((Polygon-1)*3+2):
    loadPolygonPlasma2D(D,LL,s,iteration); 
    break;
  case ((Polygon-1)*3+3):
    loadPolygonPlasma3D(D,LL,s,iteration); 
    break;
  default:
    ;
  }
}

/*
void loadMovingPlasma2DBoost(Domain *D)
{
   int i,j,istart,iend,jstart,jend,s,l,intNum,cnt,np,nc,leftIndex,rightIndex;
   int posX,posY,iter,t;
   float tmp,dx,dy,cx,cy;
   float wp,pDt,v1,v2,v3,gamma,beta[2];

   float ne,randTest=0,positionX,positionY;
   float maxwellianVelocity();
   float randomValue();

   LoadList *LL;
   ptclList *New,*p;   
   Particle **particle;
   particle=D->particle;

   dx=D->dx;
   dy=D->dy;

   beta[0]=D->beta;
   beta[1]=1.0;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   //position define      
   iter=0;
 for(i=iend-1; i<=iend; i++)
 {
   for(j=jstart; j<jend; j++)
   {
     LL=D->loadList;
     s=0;
     while(LL->next)
     {
       cx=LL->cx;
       cy=LL->cy;
       for(l=0; l<LL->lnodes-1; l++)
         for(t=0; t<LL->tnodes-1; t++)
         {
           
           if(LL->type==Circle) {
             tmp=sqrt((cx-(i-istart+D->minXSub)*dx)*(cx-(i-istart+D->minXSub)*dx)+
                    (cy-(j-jstart+D->minYSub)*dy)*(cy-(j-jstart+D->minYSub)*dy));
             posX=(int)((cx-tmp)/dx+istart)-istart;
             posY=j+D->minYSub-jstart;
           }
           else if(LL->type==Point4) {
             posX=i+D->minXSub-istart;
             posY=j+D->minYSub-jstart;
           }
 
           if((posX>=LL->lpoint[l] && posX<LL->lpoint[l+1]) && 
              (posY>=LL->tpoint[t] && posY<LL->tpoint[t+1]))
           {
             ne=((LL->ln[l+1]-LL->ln[l])/(LL->lpoint[l+1]-LL->lpoint[l])
                *(posX-LL->lpoint[l])+LL->ln[l]);
             ne*=((LL->tn[t+1]-LL->tn[t])/(LL->tpoint[t+1]-LL->tpoint[t])
                *(posY-LL->tpoint[t])+LL->tn[t]);
             ne*=LL->numberInCell*beta[iter];	//it is the float number of superparticles.
             intNum=(int)ne;
             randTest=ne-intNum;
             cnt=0;
             while(cnt<intNum)
             {               
               positionX=randomValue(beta[iter]);
               positionY=randomValue(1.0);
  
               if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
               {
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;
                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            
               } 
               else if(LL->withNextSpcs==1 && LL->withPrevSpcs==0)
               {
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;

                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            

                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s+1]->pt;
                 particle[i][j].head[s+1]->pt = New;

                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            
               } 
             
               cnt++;
             }		//end of while(cnt)

             if(randTest>randomValue(1.0))
             {
               positionX=randomValue(beta[iter]);
               positionY=randomValue(1.0);

               if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
               {
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;
 
                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            
               } 
               else if(LL->withNextSpcs==1 && LL->withPrevSpcs==0)
               {
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;
 
                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            
 
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s+1]->pt;
                 particle[i][j].head[s+1]->pt = New;
 
                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            
               } 
             }		//end of if(randTest)
           }		//end of if(l,t nodes)
         }		//end of for(lnodes,tnodes)  

       LL=LL->next;
       s++;
     }				//End of while(LL)   
   } 				//End of for(j)
   iter=iter+1;
 }				//End of for(i)
   //other define define      
   LL=D->loadList;
   s=0;
   while(LL->next)
   {
     for(i=iend-1; i<=iend; i++)
       for(j=jstart; j<jend; j++)
       {
         p=particle[i][j].head[s]->pt;
         while(p)
         {
           p->E1=p->B1=p->Pr=p->Pl=p->Sr=p->Sl=0.0;
           v1=maxwellianVelocity(LL->temperature)/velocityC;
           v2=maxwellianVelocity(LL->temperature)/velocityC;
           v3=maxwellianVelocity(LL->temperature)/velocityC;
           gamma=1.0/sqrt(1.0-(v1*v1+v2*v2+v3*v3));
           p->p1=-D->gamma*D->beta+v1;
           p->p2=v2;
           p->p3=v3;

           p=p->next;
         } 
       }				//end of for(j)
       LL=LL->next;
       s++;
   }					//End of while(LL)

}
*/

void loadPolygonPlasma2D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,t;
   float posX,posY,posZ,v1,v2,v3,centerY,centerZ,coef;
   float ne,randTest,positionX,positionY;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   Particle ***particle;
   particle=D->particle;

   ptclList *New,*p;   

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;
   k=0;
   centerY=LL->centerY;
   centerZ=LL->centerZ;
   coef=LL->coefChannel;

   srand(iteration*(myrank+1));

   //position define      
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       for(l=0; l<LL->xnodes-1; l++)
         for(t=0; t<LL->ynodes-1; t++)
         {
           posX=(float)(i+D->minXSub-istart);
           posY=(float)(j+D->minYSub-jstart);
           posZ=(float)(k+D->minZSub-kstart);
           if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1] &&
              posY>=LL->ypoint[t] && posY<LL->ypoint[t+1])
           {
             ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
             ne*=((LL->yn[t+1]-LL->yn[t])/(LL->ypoint[t+1]-LL->ypoint[t])*(posY-LL->ypoint[t])+LL->yn[t]);
             ne*=channel(centerY,posY,centerZ,posZ,coef);
             ne*=LL->numberInCell;	//it is the float number of superparticles.
             intNum=(int)ne;
             randTest=ne-intNum;
             cnt=0;
             while(cnt<intNum)
             {      
               positionX=randomValue(1.0);
               positionY=randomValue(1.0);

               New = (ptclList *)malloc(sizeof(ptclList)); 
               New->next = particle[i][j][k].head[s]->pt;
               particle[i][j][k].head[s]->pt = New;
 
               New->x = positionX;
               New->oldX=i+positionX;
               New->y = positionY;
               New->oldY=j+positionY;
               New->z = 0;
               New->oldZ=k +0;

               New->E1=New->E2=New->E3=0.0;
               New->B1=New->B2=New->B3=0.0;
               v1=maxwellianVelocity(LL->temperature)/velocityC;
               v2=maxwellianVelocity(LL->temperature)/velocityC;
               v3=maxwellianVelocity(LL->temperature)/velocityC;
               New->p1=-D->gamma*D->beta+v1;
               New->p2=v2;
               New->p3=v3;
               LL->index+=1;
               New->index=LL->index;            
               New->core=myrank;            

               cnt++; 
             }		//end of while(cnt)
/*
             if(randTest>randomValue(1.0))
             {
               positionX=randomValue(1.0);
               positionY=randomValue(1.0);

               New = (ptclList *)malloc(sizeof(ptclList)); 
               New->next = particle[i][j][k].head[s]->pt;
               particle[i][j][k].head[s]->pt = New;
 
               New->x = positionX;
               New->oldX=i+positionX;
               New->y = positionY;
               New->oldY=j+positionY;
               New->z = 0;
               New->oldZ=k +0;
               New->E1=New->E2=New->E3=0.0;
               New->B1=New->B2=New->B3=0.0;
               v1=maxwellianVelocity(LL->temperature)/velocityC;
               v2=maxwellianVelocity(LL->temperature)/velocityC;
               v3=maxwellianVelocity(LL->temperature)/velocityC;
               New->p1=-D->gamma*D->beta+v1;
               New->p2=v2;
               New->p3=v3;
               LL->index+=1;
               New->index=LL->index;            
               New->core=myrank;            
             }		//end of if(randTest)
*/
           }	
         } 		//end of for(lnodes)  
     }			//End of for(i,j)
         
}

void loadPolygonPlasma3D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,m,n;
   float posX,posY,posZ,v1,v2,v3,centerY,centerZ,coef;
   float ne,randTest,positionX,positionY,positionZ;
   Particle ***particle;
   particle=D->particle;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New,*p;   

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;
   centerY=LL->centerY;
   centerZ=LL->centerZ;
   coef=LL->coefChannel;

   srand(iteration*(myrank+1));

   //position define      
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
       for(k=kstart; k<kend; k++)
       {
         for(l=0; l<LL->xnodes-1; l++)
           for(m=0; m<LL->ynodes-1; m++)
             for(n=0; n<LL->znodes-1; n++)
             {
               posX=(float)(i+D->minXSub-istart);
               posY=(float)(j+D->minYSub-jstart);
               posZ=(float)(k+D->minZSub-kstart);
 
               if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1] &&
                  posY>=LL->ypoint[m] && posY<LL->ypoint[m+1] &&
                  posZ>=LL->zpoint[n] && posZ<LL->zpoint[n+1])
               {
                 ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
                 ne*=((LL->yn[m+1]-LL->yn[m])/(LL->ypoint[m+1]-LL->ypoint[m])*(posY-LL->ypoint[m])+LL->yn[m]);
                 ne*=((LL->zn[n+1]-LL->zn[n])/(LL->zpoint[n+1]-LL->zpoint[n])*(posZ-LL->zpoint[n])+LL->zn[n]);
                 ne*=channel(centerY,posY,centerZ,posZ,coef);
                 ne*=LL->numberInCell;	//it is the float number of superparticles.
                 intNum=(int)ne;
                 randTest=ne-intNum;
             
                 cnt=0;
                 while(cnt<intNum)
                 {               
                   positionX=randomValue(1.0);
                   positionY=randomValue(1.0);
                   positionZ=randomValue(1.0);
  
                   New = (ptclList *)malloc(sizeof(ptclList)); 
                   New->next = particle[i][j][k].head[s]->pt;
                   particle[i][j][k].head[s]->pt = New;
 
                   New->x = positionX;
                   New->oldX=i+positionX;
                   New->y = positionY;
                   New->oldY=j+positionY;
                   New->z = positionZ;
                   New->oldZ=k +positionZ;
  
                   New->E1=New->E2=New->E3=0.0;
                   New->B1=New->B2=New->B3=0.0;
                   v1=maxwellianVelocity(LL->temperature)/velocityC;
                   v2=maxwellianVelocity(LL->temperature)/velocityC;
                   v3=maxwellianVelocity(LL->temperature)/velocityC;
                   New->p1=-D->gamma*D->beta+v1;
                   New->p2=v2;
                   New->p3=v3;
                   LL->index+=1;
                   New->index=LL->index;            
                   New->core=myrank;            
   
                   cnt++; 
                 }		//end of while(cnt)

                 if(randTest>randomValue(1.0))
                 {
                   positionX=randomValue(1.0);
                   positionY=randomValue(1.0);
                   positionZ=randomValue(1.0);

                   New = (ptclList *)malloc(sizeof(ptclList)); 
                   New->next = particle[i][j][k].head[s]->pt;
                   particle[i][j][k].head[s]->pt = New;
   
                   New->x = positionX;
                   New->oldX=i+positionX;
                   New->y = positionY;
                   New->oldY=j+positionY;
                   New->z = positionZ;
                   New->oldZ=k + positionZ;
                   New->E1=New->E2=New->E3=0.0;
                   New->B1=New->B2=New->B3=0.0;
                   v1=maxwellianVelocity(LL->temperature)/velocityC;
                   v2=maxwellianVelocity(LL->temperature)/velocityC;
                   v3=maxwellianVelocity(LL->temperature)/velocityC;
                   New->p1=-D->gamma*D->beta+v1;
                   New->p2=v2;
                   New->p3=v3;
                   LL->index+=1;
                   New->index=LL->index;            
                   New->core=myrank;            
                 }		//end of if(randTest)
               }		//End of if (l,m,n)
             }			//End of for(l,m,n)

       }		//End of for(i,j,k)         
}

float maxwellianVelocity(float temperature)
{
   float vth,r,prob,v,random;
   int intRand,randRange=1e5;

   vth=sqrt(2.0*eCharge*temperature/eMass);
   
   r=1.0;
   prob=0.0;
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((float)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((float)intRand)/randRange;
      v = 6.0*(random-0.5);
      prob=exp(-v*v);
   }
   return vth*v;
}

float randomValue(float beta)
{
   float r;
   int intRand, randRange=100, rangeDev;

   rangeDev=(int)(randRange*(1.0-beta));
   intRand = rand() % (randRange-rangeDev);
   r = ((float)intRand)/randRange+(1.0-beta);

   return r;
}

