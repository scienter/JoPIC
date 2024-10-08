#define Electron 	1
#define Test	 	2
#define HPlus0 	 	100
#define HPlus1 	 	101
#define HePlus0 	200
#define HePlus1 	201
#define HePlus2 	202
#define CPlus0          600
#define CPlus1          601
#define CPlus2          602
#define CPlus3          603
#define CPlus4          604
#define CPlus5          605
#define CPlus6          606
#define NPlus0          700
#define NPlus1          701
#define NPlus2          702
#define NPlus3          703
#define NPlus4          704
#define NPlus5          705
#define NPlus6          706
#define NPlus7          707
#define NePlus0         1000
#define NePlus1         1001
#define NePlus2         1002
#define NePlus3         1003
#define NePlus4         1004
#define NePlus5         1005
#define NePlus6         1006
#define NePlus7         1007
#define NePlus8         1008
#define NePlus9         1009
#define NePlus10        1010
#define AlPlus4         1304
#define userDefined   	9999999

#define Polygon    	1
#define Defined    	2
#define Beam	    	3
#define Channel    	4
#define BoostFrame	5
#define Circle    	6
#define Exp    		7

#define Constant   	0
#define Gaussian   	1
#define Polynomial   	2

#define byNumber	0
#define byDensity	1

typedef struct _LoadList  {
   int type;
   int species;
   double superP;
   double density;
   double numberRZ;
   double numberPhi;
   double criticalDensity;
   double targetW;
   double num;      //exceeded number of particle which is less than 1
   int xnodes;     //longitudinal point number
   double *xn;      //longitudinal density (1 is P->density)
   double *xpoint;    //longitudinal point
   int ynodes;     //transverse point number
   double *yn;      //transverse density (1 is P->density)
   double *ypoint;    //transverse point
   int znodes;     //transverse point number
   double *zn;      //transverse density (1 is P->density)
   double *zpoint;    //transverse point
   double givenMinPx;	//for saveParticle option
   int ChXnodes;     //longitudinal point number for chaanel
   double *ChXn;      //longitudinal coef. value for channel
   double *ChXpoint;    //longitudinal point for channel

   //initial momentum distribution
   double z0,pz0;
   double delZ;

   //defined plasma
   int defineMode;
   double *xPosition;
   double *yPosition;
   double *zPosition;
   int numDefined;
   double xLengDef;
   double yLengDef;
   double zLengDef;
   int numDefPtcls;
   double **define;
   int maxLoadTime;
   int minLoadTime;
   double minX;
   double maxX;
   double minY;
   double maxY;
   double minZ;
   double maxZ;

   int pointPosition;   
   double pz,px,py;

   double mass;
   int charge;
   
   double temperature;
   
   //channel
   double channelCoef;

   //ionization
   int levels;
   int ionFinal;
   double givenMinA0;
   double *ionEnergy;
   double *ionW;

   // Beam
	int loadingStep,gaussMode,numGam;
	double energy,peakCurr,eChirp,posZ,spread;
	double gaussPower,sigZ,emitR,betaR,alphaR,focalL;

   //pair current
   int pair;

   struct _LoadList *next;
} LoadList;

typedef struct _PlasmaLens  {
	int xnodes;       //longitudinal point number
	double *xn;       //longitudinal density (1 is P->density)
	double *xpoint;      //longitudinal point
	double radius;
	double current;

	struct _PlasmaLens *next;
} PlasmaLens;

