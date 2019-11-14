# define _BSD_SOURCE
# define _XOPEN_SOURCE 
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <math.h>
# include <string.h>
# include <unistd.h>
# include <limits.h>
# include <sys/types.h>
# include <stdint.h>
# include <ctype.h>
# include <hdf5.h>
# include <hdf5_hl.h>

void date_from_sec(char *date,time_t secs);
time_t sec_from_date(char *date);
int H5get_variable_string(hid_t h5,char *datasetname, char *attrname, char *str);
uint8_t dBZIfromAcc(int32_t Acc);
void gen_dBZIfromIR(void);

double ZR_A,ZR_B,AccScaler,gain,offset;
uint16_t nodata,newNodata;
/* uint8_t *dBZIfromIR; */
uint8_t dBZIfromIR[65536]={0};
uint8_t GENERATE_DBZ=0, GENERATE_INTER_ACC=0;

int main(int argc, char *argv[])
{

   /* HDF5 variables */
  hid_t INH5;
  hid_t datatype;
  hsize_t dims[3];
  H5T_class_t class_id;

  size_t typesize;
  long SAMPLE = 1044550;
  uint8_t *PGMAccdata,**dBZIdata;
  uint16_t obsZItoncZI[65536]={0};
  uint16_t *ZIarr0,*ZIarr1; 
  int32_t ***AccAverLUT,*Accdata,PGMAcc,PGMscaler=1000,PGMscaler2,intlen;
  int32_t *AccFromZI,ZIveclen,**trajecLUT[2];
  int i,timesteps=0,I,intsteps=10,nodata;
  char *p=NULL,*timestamp,datasetname[200],pidname[255];
  char *membername,*h5file,*obsfile,*outdir,varstr[1000];
  long xsize,ysize,csize,arrsize,fieldsize;
  long N;
  double *motarr,vecscaler;
  double fnodata;
  time_t obsecs,dT,IdT;
  pid_t pid;
  FILE *PIDF;

  pid=getpid();
  setbuf(stdout,NULL);
  setbuf(stderr,NULL);
  H5Eset_auto(H5E_DEFAULT,NULL,NULL);

  timestamp=argv[1];
  membername=argv[2];
  h5file=argv[3];
  obsfile=argv[4];
  outdir=argv[5];
  if(argv[6]) intsteps=atoi(argv[6]);
  vecscaler=1.0/(double)intsteps;
  intlen=intsteps-1;

  if(p=getenv("PGMSCALER")) PGMscaler=atol(p);
  if(p=getenv("TIMESTEPS")) timesteps=atoi(p);
  if(p=getenv("SAMPLE_INDEX")) SAMPLE=atol(p); else SAMPLE=0;
  if(getenv("GENERATE_DBZ")) GENERATE_DBZ=1;  
  if(getenv("GENERATE_INTER_ACC")) GENERATE_INTER_ACC=1;  

  sprintf(pidname,"%s_%s.pid",timestamp,membername);
  PIDF = fopen(pidname,"w");
  fprintf(PIDF,"%d",(int)pid);
  fclose(PIDF);

  PGMscaler2=2*PGMscaler;
  obsecs=sec_from_date(timestamp);

  /* Open nowcast file and read configuration attributes */
  INH5=H5Fopen(h5file, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(H5get_variable_string(INH5,"/meta/configuration","ZR_A",varstr) >= 0 )
     ZR_A=atof(varstr); else ZR_A=223.0;
  if(H5get_variable_string(INH5,"/meta/configuration","ZR_B",varstr) >= 0 )
     ZR_B=atof(varstr); else ZR_B=1.53;
  if(H5get_variable_string(INH5,"/meta/configuration","NOWCAST_TIMESTEP",varstr) >= 0 )
     dT=60*atoi(varstr); else dT=300;
  if(H5get_variable_string(INH5,"/meta/configuration","NUM_TIMESTEPS",varstr) >= 0 )
    if(!timesteps) timesteps=atoi(varstr); /* 6h with 5min steps */

  IdT=dT/intsteps; /* interpolation timestep */

  /* Read data conversion attributes */  
  sprintf(datasetname,"/member-00/leadtime-00");
  H5LTget_attribute_double(INH5,datasetname,"gain",&gain);
  H5LTget_attribute_double(INH5,datasetname,"nodata",&fnodata);
  H5LTget_attribute_double(INH5,datasetname,"offset",&offset);
  nodata=(uint16_t)(fnodata+1e-6);

  /* printf("%f %f %f %f %d\n",ZR_A, ZR_B, gain, offset, nodata); */

  /* Read the motion field */
  sprintf(datasetname,"/%s/motion",membername);
  /* sprintf(datasetname,"/motion"); */
  H5LTget_dataset_info(INH5,datasetname, dims, &class_id, &typesize );
  csize=(long)dims[0];
  ysize=(long)dims[1];
  xsize=(long)dims[2];
  fieldsize=xsize*ysize;
  arrsize=fieldsize*csize;
  motarr=malloc(arrsize*typesize);
  ZIarr0=malloc(fieldsize*sizeof(uint16_t));
  ZIarr1=malloc(fieldsize*sizeof(uint16_t));
  datatype=H5T_IEEE_F64LE; /* for motion data */
  H5LTread_dataset(INH5,datasetname,datatype,motarr);
  gen_dBZIfromIR();
  for(N=0;N<arrsize;N++) motarr[N]*=vecscaler;

  /* Create LUT for accumulation per one nowcast timestep for each dBZI value pair */
  {
     double B,C,R,dBZ,Rscaler,k;
     int32_t zN,zN0,zN1,maxdBZ=100,maxN; /* maxN = 1320, jos gain=0.1 ja offset=-32 */
     int32_t Acc,Acc0,Acc1;
     int32_t X,rX,X1;
     int32_t IR; /* [1/100000 mm] */
  

     maxN = (maxdBZ-offset)/gain;
     
     newNodata = maxN+1;
     ZIveclen = newNodata+1;

     AccFromZI = calloc(ZIveclen,sizeof(int32_t));
     AccAverLUT = malloc(ZIveclen*sizeof(int32_t **));
    /* Conversion from intensity [mm/h] to accumulation during one IdT (interp timestep) 
       unit 1/100000 mm/IdT */
     Rscaler = dT/0.036/(double)intsteps; 
     AccScaler = Rscaler/269.0; 
     B = 0.1/ZR_B;
     C = log10(ZR_A)/ZR_B;

     /* zN -> Acc(dT) vector and memory allocations */

     /* printf("%d %f %d %d\n",intsteps, Rscaler, maxN, ZIveclen); */

     for(zN0=0;zN0<ZIveclen;zN0++)
     {
        AccAverLUT[zN0] = malloc(ZIveclen*sizeof(int32_t *));
        for(zN1=0;zN1<ZIveclen;zN1++) AccAverLUT[zN0][zN1]=calloc(intlen,sizeof(int32_t));        
        if(zN0)
	{
           dBZ = gain*(double)zN0+offset;
           R = pow(10.0,B*dBZ - C);
           IR = (int32_t)(R*Rscaler);
           AccFromZI[zN0]=IR;
	}
     }
     AccFromZI[maxN]=-1;
     /* printf("TOIMII 1\n"); */
  
     /* Either value nodata */
     for(zN=0;zN<ZIveclen;zN++) for(X=0 ; X<intlen ; X++) 
         AccAverLUT[newNodata][zN][X] = AccAverLUT[zN][newNodata][X] = -1; 

     /* Linear interpolation of two timestep-accumulation over interp steps */
     for(zN0=1; zN0<ZIveclen; zN0++) 
     {
        Acc0=AccFromZI[zN0];
        for(zN1=zN0; zN1<ZIveclen; zN1++) 
        {
           Acc1=AccFromZI[zN1];
           k=(Acc1-Acc0)/(double)intsteps;
           for(X=0,X1=1,rX=intlen-1 ; X<intlen ; X++,X1++,rX--)        
	   {  
              if(zN0 == zN1) AccAverLUT[zN0][zN0][X] = Acc0;
              else
	      { 
		 Acc = (int32_t)(Acc0 + k*(double)X1);
                 AccAverLUT[zN0][zN1][X] = AccAverLUT[zN1][zN0][rX] = Acc;
              }
	   }
	}
     }
  }

  /* Create obsdata to ncdata LUT */
  {
     uint16_t N;
 
     /*
     obsZItoncZI[255] = newNodata;
     obsZItoncZI[65535] = newNodata;
     */
     obsZItoncZI[255] = 0;
     obsZItoncZI[65535] = 0;
     for(N=1;N<255;N++) obsZItoncZI[N] = (uint16_t)((0.5*(double)N-32.0-offset)/gain); /* 1-byte IRIS dBZ */
     for(N=256;N<65535;N++) obsZItoncZI[N] = (uint16_t)((0.01*(double)N-327.68-offset)/gain); /* 2-byte IRIS dBZ */
  }


  /* printf("TOIMII 2\n"); */
  

  /* Create N-LUTs for "from past" and "from future" motion for each interpolation timestep index */
  {
     double dX[2],dY[2];
     long N,motN,X,Y,iX,iY,iI,iN;
     short T,sig,PAST=0,FUTU=1,wI[2];
     
  
     /* Source point trajectory arrays */
     for(T=PAST;T<=FUTU;T++)
     {
        trajecLUT[T]=calloc(intlen,sizeof(int32_t *));
        for(iI=0 ; iI<intlen ; iI++)
        {
           trajecLUT[T][iI]=malloc(fieldsize*sizeof(int32_t));
        }
     }
     /* printf("TOIMII 3a\n"); */

     /* Loop thru destination area */
     for(Y=0;Y<ysize;Y++) for(X=0;X<xsize;X++)
     {
        N=Y*xsize+X; 
 
        /* First step of trajectories (1/intstep of original motion vector) */
        dX[PAST] = dX[FUTU] = motarr[N];
        dY[PAST] = dY[FUTU] = motarr[N + fieldsize];

        /* Interpolation loop, begins one step from past field (obs field at start) */
        for(iI=0 ; iI<intlen ; iI++)
        {
	  /* weight (distance) indices */
	   wI[PAST]=iI;
           wI[FUTU]=intlen-iI-1;

	   for(T=PAST,sig=-1; T<=FUTU ; T++,sig=1) /* trajectory "from past" has reverse motion field */
	   {
	      iX = (long)((double)X+0.5 + (double)sig*dX[T]);
	      iY = (long)((double)Y+0.5 + (double)sig*dY[T]);

	      /* If source pixel outside area, set trajectory source point array index to -1 */
              if(iX>=xsize ||
                 iX<0      ||
                 iY>=ysize ||
                 iY<0        ) { iN = -1; motN = N; } 
              else iN = motN = iY*xsize + iX;
	      /*
              if(SAMPLE) if(N==SAMPLE)
	      {
		printf("%d\t%d\t\tdX=%f\tdY=%f\tiX=%ld\tiY=%ld\tiN=%ld\n",T,iI,dX[T],dY[T],iX,iY,iN); 
	      } 
	      */             
              /* trajectory increment */
              dX[T] += motarr[motN];
              dY[T] += motarr[motN + fieldsize];
              /* Assign source point index iN and point's distance (weight) wI to destination LUT */
              trajecLUT[T][wI[T]][N] = iN;
	   }
           /* if(SAMPLE) if(N==SAMPLE) printf("\n"); */
	}
     }
  }


  /* printf("TOIMII 3\n"); */

  /* Read last observation dBZ data */
  {
     FILE *INPGM;
     int obsdyn=255,N;
     uint8_t *obsdata8;
     uint16_t dBZI,OBS8=0,*obsdata16;
     char hdr[256];     

     INPGM=fopen(obsfile,"r");
     for(i=0;i<3;i++)
     {
        memset(hdr,0,256);
        fgets(hdr,255,INPGM);
        if(hdr[0]=='#') { i--; continue; }
        if(i==2) obsdyn=atoi(hdr);
    /*    if(i==1) sscanf(hdr,"%ld %ld",&Xdim,&Ydim); */
     }
     /* printf("TOIMII 4\n"); */

     if(obsdyn>255)
     { 
        obsdata16=malloc(fieldsize*2);
        fread(obsdata16,1,fieldsize*2,INPGM);
        swab(obsdata16,obsdata16,fieldsize*2);
        OBS8=0;
     } 
     else
     { 
        obsdata8=malloc(fieldsize);
        fread(obsdata8,1,fieldsize,INPGM);
        OBS8=1;
     } 
     fclose(INPGM);

     /* printf("TOIMII 5\n"); */

    /* Convert obsdata to nc-type data */
     for(N=0;N<fieldsize;N++)
     { 
        if(OBS8) dBZI=(uint16_t)obsdata8[N]; else dBZI=obsdata16[N];        
        ZIarr0[N] = obsZItoncZI[dBZI];
     }
    
     if(OBS8) free(obsdata8); else free(obsdata16);
  }

  PGMAccdata = malloc(fieldsize);
  Accdata = calloc(fieldsize,sizeof(int32_t));
  dBZIdata=malloc(intsteps*sizeof(uint8_t *));
  for(I=0;I<intsteps;I++) dBZIdata[I]=malloc(fieldsize);

  /* printf("TOIMII 6\n"); */




/* Main time loop ________________________________________________________________________________________________ */

  {
     int tI,iI,wI;
     int32_t pastAcc,Acc; /* futuAcc; */
     uint16_t pastZI,futuZI,*tmpPtr,*pastZIarrPtr,*futuZIarrPtr; 
     long N,pastN,futuN,X,Y;
     time_t ncsecs;
     char nowcstamp[15],intstamp[15];
     char outpath[250], datasetname[200];
     hid_t datatype=H5T_STD_U16LE; /* For input nc data */
     FILE *OUTFILE;

     /* Assign pointers to nowcast data arrays (integer dBZ values) */
     pastZIarrPtr = ZIarr0; /* in the beginning, this points to obs data */
     futuZIarrPtr = ZIarr1;

     for(tI=0;tI<timesteps;tI++)
     {
        ncsecs=dT*(tI+1);
        date_from_sec(nowcstamp,obsecs+ncsecs);

        /* Read nc data for timestep tI to future data array */
        sprintf(datasetname,"/%s/leadtime-%02d",membername,tI);
        H5LTread_dataset(INH5,datasetname,datatype,futuZIarrPtr);
        for(N=0;N<fieldsize;N++) if(futuZIarrPtr[N] == nodata) futuZIarrPtr[N]=newNodata;

        printf("Data tI=%d of %s read\n",tI,datasetname);

        /* reset monitor accumulation field */
        memset(PGMAccdata,0,fieldsize);
        for(iI=0;iI<intsteps;iI++) memset(dBZIdata[iI],0,fieldsize);

        /* printf("INIT tI %d OK\n",tI); */

        /* Loop thru area */
        for(Y=0;Y<ysize;Y++) for(X=0;X<xsize;X++)
        {
	   N = pastN = Y*xsize+X; /* pixel index of "real time" data */

           /* Interpolation loop, begins witth past data */ 
           for(iI=0,wI=-1; iI<intsteps ; iI++,wI++)        
	   {
	      /* data pixel from past data */ 
	      if(wI>=0) pastN = trajecLUT[0][wI][N];   /* trajectory starting point from past data  */

              if(pastN>=0)
	      {
                 pastZI = pastZIarrPtr[pastN]; /* integer dBZ */
                 pastAcc = AccFromZI[pastZI];  /* conversion from dBZ to accumulation for one interpolation time step */

                 if(iI)
		 {
	            /* data pixel from future data */ 
                    futuN = trajecLUT[1][wI][N];
                    if(futuN>=0)
		    {
                       futuZI = futuZIarrPtr[futuN];
                       /* futuAcc = AccFromZI[futuZI]; no need to fetch */

                       /* time-weighted mean of past and future accumulation */
	               Acc  = AccAverLUT[pastZI][futuZI][wI];
		    } else Acc = -1;
		 } else Acc = pastAcc;
	      } else Acc = -1;

              if(GENERATE_DBZ) dBZIdata[iI][N] = dBZIfromAcc(Acc);
              Accdata[N]+=Acc;
# if 0
              if(SAMPLE) if(N==SAMPLE)
	      { 
		/* AAcc = pastAcc+iI*(futuAcc-pastAcc)/intsteps; */
                  printf("N: %ld %ld %ld\n",N,pastN,futuN);
		/*  printf("A: %ld\tpastZI=%04ld\tfutuZI=%04ld\tpastAcc=%04d\tfutuAcc=%04d\t\tAcc=%d\n",iI, pastZI, futuZI, pastAcc, futuAcc, Acc); */
                 /* printf("M: %ld\ttX=%04ld\ttY=%04ld\tsX=%04ld\tsY=%04ld\t\t%.3f\t%.3f\n",iI,tX,tY,sX,sY,fu,fv); */
	      }
# endif
           }
        }

        /* swap nowcast data array pointers to move future data to past data for next time step */
        tmpPtr = pastZIarrPtr;
        pastZIarrPtr = futuZIarrPtr;
        futuZIarrPtr = tmpPtr;
     
        printf("tI=%d %s interpolated\n",tI,nowcstamp);        

        /* create PGM data for monitoring if wanted */
        if(GENERATE_DBZ) for(iI=0;iI<intsteps;iI++)
        {
           date_from_sec(intstamp,obsecs + tI*dT + iI*IdT);
           sprintf(outpath,"%s/Z_%s_%s-%s.pgm",outdir,membername,timestamp,intstamp);
           OUTFILE=fopen(outpath,"w");
           fprintf(OUTFILE,"P5\n%ld %ld\n255\n",xsize,ysize);
           fwrite(dBZIdata[iI],1,fieldsize,OUTFILE);
           fclose(OUTFILE);
           printf("dBZ data %s written\n",intstamp);
        }

        if(GENERATE_INTER_ACC) 
	{
           for(N=0;N<fieldsize;N++)
           { 
              if(Accdata[N]<0) { PGMAccdata[N]=255; continue; }
              if(! Accdata[N]) continue;
              if(Accdata[N]<PGMscaler2) { PGMAccdata[N]=1; continue; }
              PGMAcc = Accdata[N]/PGMscaler;
              if(PGMAcc>250) PGMAccdata[N]=250; else PGMAccdata[N]=PGMAcc;
           }

           sprintf(outpath,"%s/Interacc_%s_%s-%s.pgm",outdir,membername,timestamp,nowcstamp);
           OUTFILE=fopen(outpath,"w");
           fprintf(OUTFILE,"P5\n%ld %ld\n255\n",xsize,ysize);
           fwrite(PGMAccdata,1,fieldsize,OUTFILE);
           fclose(OUTFILE);
	}
     

        /* Write RAVAKE-style accumulation data for one member only */
     
        sprintf(outpath,"%s/%s_RAVACC_%s-%.12s+%03d.dat",
                outdir,membername,timestamp,nowcstamp,(int)(ncsecs/60));
        OUTFILE=fopen(outpath,"w");
        fwrite(Accdata,1,fieldsize*sizeof(int32_t),OUTFILE);
        fclose(OUTFILE);
        printf("Nowcast data %s processed\n",nowcstamp);
     }
  }

/* ________________________________________________________________________________________________ */

  {
     int i,j;

     for(i=0;i<ZIveclen;i++)
     { 
        for(j=0;j<ZIveclen;j++) free(AccAverLUT[i][j]); 
        free(AccAverLUT[i]);
     }

     for(i=0;i<intsteps;i++)
     { 
        free(dBZIdata[i]);
        if(i<intlen)
	{
           free(trajecLUT[0][i]);
           free(trajecLUT[1][i]);
	}
     }
     free(dBZIdata);
     free(trajecLUT[0]);
     free(trajecLUT[1]);
     free(AccAverLUT);
     /* printf("TOIMII 7\n"); */
     free(AccFromZI);
     free(PGMAccdata);
     free(Accdata);
     free(ZIarr0);
     free(ZIarr1);
     free(motarr);
  }

  remove(pidname);
  return(0);
}

/* ################################################################################################ */

void date_from_sec(char *date,time_t secs)
{
   struct tm *Sdd;

   Sdd=gmtime(&secs);
   sprintf(date,"%d%02d%02d%02d%02d%02d",Sdd->tm_year+1900,Sdd->tm_mon+1,
                                         Sdd->tm_mday,Sdd->tm_hour,
	                                 Sdd->tm_min,Sdd->tm_sec);
   return;
}

/* ================================================================================================ */

time_t sec_from_date(char *date)
{
   struct tm Sd,*Sdd;
   int y,m;
   time_t secs;
 
   sscanf(date,"%4d%2d%2d%2d%2d",&y,&m,&Sd.tm_mday,&Sd.tm_hour,&Sd.tm_min);
   Sd.tm_year=y-1900;
   Sd.tm_mon=m-1;
   Sd.tm_isdst=0;
   Sd.tm_sec=0;
   Sdd=&Sd;
   secs=mktime(Sdd);
   return(secs);
}

/* ================================================================================================ */

int H5get_variable_string(hid_t h5,char *datasetname, char *attrname, char *str)
{
   hid_t dataset,attr,memtype,space;
   hvl_t  rdata;             /* Pointer to vlen structures */
   char *ptr;
   int i,len;
   /*
   int ndims,i,len;
   hsize_t dims[1] = {1};
   */
   herr_t status;

   memset(str,0,strlen(str));
   dataset=H5Dopen(h5,datasetname,H5P_DEFAULT);
   if(dataset<0) dataset=H5Gopen(h5,datasetname,H5P_DEFAULT);
   if(dataset<0) return(-1);
   attr = H5Aopen(dataset, attrname, H5P_DEFAULT);

   space = H5Aget_space(attr);
   memtype = H5Tvlen_create(H5T_NATIVE_CHAR);
   status = H5Aread(attr, memtype, &rdata);
   if(status<0) return(-1);
   ptr = rdata.p;
   len = rdata.len;
   for (i=0; i<len; i++) str[i]=ptr[i];

   status = H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, &rdata);
   status = H5Aclose (attr);
   status = H5Dclose (dataset);
   status = H5Sclose (space);
   status = H5Tclose (memtype);
   return(len);
}

uint8_t iRtodBZI(int32_t IR)
{
   int ZI;
   uint8_t dBZI;
   double R,dBZ; 

   if(!IR) dBZI=0;
   else
   {
      R=(double)IR;
      dBZ=10.0*log10(ZR_A*pow(R,ZR_B));
      if(dBZ > -32.0) ZI=(int)(2.0*dBZ)+64; else ZI=0;
      if(ZI > 254) ZI=254;
      dBZI=(uint8_t)ZI;
   }

    /*    printf("%d=%d ",R,dBZI); */
   return(dBZI);
}   

uint8_t dBZIfromAcc(int32_t Acc)
{
   int32_t IR;

   if(Acc<0) return(255);
   IR=(int32_t)((double)Acc/AccScaler);
   if(!IR) return(0);
   if(IR>=65535) return(185);
   return(dBZIfromIR[IR]);
}

void gen_dBZIfromIR(void)
{
   /* construct scaled R -> dBZI LUT here */
   int ZI;
   int32_t IR;
   uint8_t dBZI;
   double R,dBZ; 

   for(IR=1;IR<65535;IR++)
   {
      R=(double)IR/269.0;
      dBZ=10.0*log10(ZR_A*pow(R,ZR_B));
      if(dBZ > -32.0) ZI=(int)(2.0*dBZ)+64; else ZI=0;
      dBZI=(uint8_t)ZI;
      dBZIfromIR[IR]=dBZI;
   }
}   
