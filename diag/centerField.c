#include "stdio.h"
#include "stdlib.h"
#include "math.h"

int main(int argc, char *argv[])
{
   double x,y,Ex,Ey,Ez,Bx,By,Bz,den,D,F,dx,dt,minX,maxX,weight,rangeX,oldEx;
	double V,maxV,posX;
   double *dataX,*dataT,***dataF,**dataD,*data1,*data2,*data3;
   FILE *in,*out;
   char fileName[100],fileName1[100],fileName2[100],outFile[100],dataType[100];
   int dimension,mode,core,initial,final,saveStep;
   int cnt,cntT,i,j,k,rangeI,step,index,cenFlag;

   if(argc < 6)   { 
      printf("centerfield mode dimension initial final saveStep dataType\n");
      printf("mode 0(x-T) : rangeX dtRatio dx core\n");
      printf("mode 1(x-T,cylind) : rangeX dtRatio dx core\n");
      printf("mode 2(0density) : \n");
      exit(0);
   }

   mode = atoi(argv[1]);
   dimension = atoi(argv[2]);
   initial = atoi(argv[3]);
   final = atoi(argv[4]);
   saveStep = atoi(argv[5]);
   cntT=(final-initial)/saveStep+1;

   sprintf(dataType,"%s",argv[6]);

   switch (mode) {
   case 0 :
   case 1 :
     rangeX = atof(argv[7]);
     dt = atof(argv[8]);
     dx = atof(argv[9]);
     core = atoi(argv[10]);

     step=initial;
//     if(mode==0) {
//       if(dimension==1) sprintf(fileName,"cenfieldE%d_Ex",step); 
//       else             sprintf(fileName,"cenYee%d_%d_XY",step,core); 
//       else             sprintf(fileName,"cenPukhov%d_%d",step,core); 
//     } else {      
       sprintf(fileName,"%s%d_%d",dataType,step,core); 
//     }

     if(fopen(fileName,"r")==NULL) {
       printf("%s is not exited.\n",fileName);
       cnt=1;
       step=final+saveStep;
     } else {
       in = fopen(fileName,"r");
       cnt=0;
       if(dimension==1) { 
         while(fscanf(in,"%lf %lf",&x,&Ex)!=EOF) cnt++;
       } else {
         while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&Ex,&Ey,&Ez,&Bx,&By,&Bz,&den)!=EOF) 
           cnt++;
       }
       fclose(in);
       step=initial;
		 cnt=(int)(rangeX/dx);
     }

     minX=-rangeX;
     rangeI=(int)(rangeX/dx)+1;

     dataX=(double *)malloc(cnt*sizeof(double));
     data1=(double *)malloc(cnt*sizeof(double));
     data2=(double *)malloc(cnt*sizeof(double));
     data3=(double *)malloc(cnt*sizeof(double));
     dataT=(double *)malloc(cntT*sizeof(double));
     dataF=(double ***)malloc(cntT*sizeof(double **));
     dataD=(double **)malloc(cntT*sizeof(double *));
     for(i=0; i<cntT; i++) {
       dataF[i]=(double **)malloc(cnt*sizeof(double *));
       for(j=0; j<cnt; j++) 
         dataF[i][j]=(double *)malloc(2*sizeof(double ));
       dataD[i]=(double *)malloc(cnt*sizeof(double));
     }

     for(j=0; j<cntT; j++) 
       for(i=0; i<cnt; i++) {
         for(k=0; k<2; k++) dataF[j][i][k]=0.0;
         dataD[j][i]=0.0;
       }

     while(step<=final) {
       j=(step-initial)/saveStep;
       dataT[j]=step;

		 maxX=step*dx*dt;
		 minX=maxX-rangeX;

       if(dimension==1) {
         sprintf(fileName,"cenfieldE%d_Ex",step);
         sprintf(fileName1,"cenfieldE%d_Ey",step);
         sprintf(fileName2,"cen0density%d_0",step);
         if(fopen(fileName,"r")==NULL) printf("%s is not exited.\n",fileName);
         if(fopen(fileName1,"r")==NULL) printf("%s is not exited.\n",fileName1);
         if(fopen(fileName2,"r")==NULL) printf("%s is not exited.\n",fileName2);
       }
       else { 
//         if(mode==0)   sprintf(fileName,"cenPukhov%d_%d",step,core);
//         if(mode==0)   sprintf(fileName,"cenYee%d_%d",step,core);
//         else          sprintf(fileName,"%s%d_%d",dataType,step,core);
//         if(fopen(fileName,"r")==NULL) printf("%s is not exited.\n",fileName);
          sprintf(fileName,"%s%d_%d",dataType,step,core);
       }
  
       
       if(dimension==1) {
         in = fopen(fileName,"r");
         for(i=0; i<cnt; i++) {
           fscanf(in,"%lf %lf",&x,&data1[i]); dataX[i]=x-dx*dt*step; }
         fclose(in);

         in = fopen(fileName1,"r");
         for(i=0; i<cnt; i++) fscanf(in,"%lf %lf",&x,&data2[i]); 
         fclose(in);

         in = fopen(fileName2,"r");
         for(i=0; i<cnt; i++) fscanf(in,"%lf %lf",&x,&data3[i]); 
         fclose(in);
       } else {
         in = fopen(fileName,"r");
			while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&Ex,&Ey,&Ez,&Bx,&By,&Bz,&den)!=EOF) {
				i=(int)((x-minX)/dx);
				if(i>=0 && i<cnt) {
           		dataF[j][i][0]=Ex;
           		dataF[j][i][1]=Ey;
           		dataD[j][i]=den;
				} else ;
         }
         fclose(in);
       }
       printf("%s is done\n",fileName);
      
       step+=saveStep;
     }

     sprintf(fileName,"cenXT");
     out = fopen(fileName,"w");
     for(j=0; j<cntT; j++) {
       for(i=0; i<cnt; i++) {
         x=-rangeX+i*dx;
         fprintf(out,"%.10g %lf %lf %lf %lf\n",x,dataT[j],dataF[j][i][0],dataF[j][i][1],dataD[j][i]);
       }
       fprintf(out,"\n");
     }
     fclose(out);    
     printf("%s is made.\n",fileName);

     sprintf(fileName,"maxTrack");
     out = fopen(fileName,"w");
     for(j=0; j<cntT; j++) {
			maxV=-1e6;
			posX=0.0;
      	for(i=0; i<cnt; i++) {
				V=fabs(dataF[j][i][1]);
				if(V>=maxV) { maxV=V; posX=-rangeX+i*dx; } else ;
		 	}
			fprintf(out,"%.10g %lf\n",posX,dataT[j]);
	  }
     fclose(out);    
     printf("%s is made.\n",fileName);
			

     
     for(i=0; i<cntT; i++)  {
       for(j=0; j<rangeI; j++) free(dataF[i][j]);
       free(dataF[i]);
       free(dataD[i]);
     }
     free(dataF); free(dataD); 
     free(dataX); free(dataT);
     free(data1); free(data2); free(data3);
     break;


   case 2 :
     step=initial;
     while(step<=final) {
       sprintf(outFile,"cenDen%d",step);
       out = fopen(outFile,"w");

       sprintf(fileName,"0density%d_0_XY",step);
       in = fopen(fileName,"r");
       while(fscanf(in,"%lf %lf %lf",&x,&y,&den)!=EOF) {
         if(y==0) fprintf(out,"%d %g %g\n",step,x,den); else;
       }
       fclose(in);
       fclose(out);
       printf("%s is done\n",fileName);
       step+=saveStep;
     }
     break;
   }

   return 0;

}
