

typedef struct _ptclHead  {
    struct _ptclList *pt;
}   ptclHead;

typedef struct _ptclList  {
    double z,oldZ,oldZ2; 
    double x,oldX,oldX2;   
    double y,oldY,oldY2;   
    double pz,px,py;
    double Ez,Ex,Ey,Bz,Bx,By;    
    double weight,charge;    
    int index,core; 
    struct _ptclList *next;
} ptclList;

