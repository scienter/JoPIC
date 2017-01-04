#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include "hdf5.h"
 
// The number of cells in the X, Y dimensions
#define NX 30
#define NY 20
 
void
write_xdmf_xml(char *name,int nx,int ny)
{
    FILE *xmf = 0;
    char filename[100];
 
    printf("here is xdmf function, %s\n",name);
    //
     // Open the file and write the XML description of the mesh..
     //
    sprintf(filename,"%s.xmf",name);
    xmf = fopen(filename,"w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
    fprintf(xmf, "     <Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n", ny,nx);
    fprintf(xmf, "     <Geometry GeometryType=\"VXVY\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",nx);
    fprintf(xmf, "        %s.h5:/X\n",name);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",ny);
    fprintf(xmf, "        %s.h5:/Y\n",name);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    fprintf(xmf, "     <Attribute Name=\"Ex\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",ny,nx);
    fprintf(xmf, "        %s.h5:/Ex\n",name);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "     <Attribute Name=\"Ey\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n",ny,nx);
    fprintf(xmf, "        %s.h5:/Ey\n",name);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}
 
int
main(int argc, char *argv[])
{
    char name[100];
    int nx,ny;

    if(argc<3)
    {
      printf("filename(.h5) nx ny\n");
      exit(0);
    }
    
    sprintf(name,"%s",argv[1]);
    nx=atoi(argv[2]);
    ny=atoi(argv[3]);

    write_xdmf_xml(name,nx,ny);
 
    return 0;
}
