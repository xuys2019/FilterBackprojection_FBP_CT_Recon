#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define pi 3.14159265358979

typedef struct TestStruct {
    float    ScanR;
    float    DistD;
    int       YL;
    int       ZL;
    float    dectorYoffset;
    float    dectorZoffset;
    float    XOffSet;
    float    YOffSet;
    float    ZOffSet;
    float    phantomXOffSet;
    float    phantomYOffSet;
    float    phantomZOffSet;
    float    DecFanAng;
    float    DecHeight;
    float    DecWidth;
    float    h;
    float    BetaS;
    float    BetaE;
    int       AngleNumber;
    int       N_2pi;
    float    Radius;
    int       RecSize;
    int       RecSizeZ;
//    float    sliceThickness;
//    int       FOILength;
//    int       FOIWidth;
//    int       FOIHeight;
    float    delta;
    float    HSCoef;
    float    k1;
    float    ***GF;
    float    ***RecIm;
} TestStruct;

extern void fbp(TestStruct *t) {

	float ScanR, DistD,DecL,DecHeight,DecWidth,ObjR,sliceThickness,dectorYoffset,dectorZoffset;
	float dx,dy,dz,dYL,dZL,DeltaFai,YLC,ZLC,XOffSet,YOffSet,ZOffSet;
	float XNC, YNC, ZNC, phantomXOffSet,phantomYOffSet,phantomZOffSet;
	float h,h1, BetaE, BetaS, delta, HSCoef, k1;
	int YL,ZL,PN,RecSize, RecSizeZ, N_2pi,N_pi,XN, YN, ZN;
//	int FOILength,FOIWidth,FOIHeight;

	ScanR = t->ScanR;        /*source object distance*/
    DistD = t->DistD;
	DecL = t->DecFanAng;     /*Project Angle*/
	YL = t->YL;              /*Detection Number */
    ZL = t->ZL;
    DecHeight = t->DecHeight;
    DecWidth = t->DecWidth;
    h = t->h;
    h = h*DecHeight;
	ObjR = t->Radius;      /*Radius of the field of view*/
	RecSize = t->RecSize;  /*reconstruction size x*/
	RecSizeZ = t->RecSizeZ;  /*reconstruction size x*/
	delta = t->delta;
	HSCoef = t->HSCoef;
	k1 = t->k1;

    BetaS = t->BetaS;
    N_2pi = t->N_2pi;
    PN = t->AngleNumber;/*Projection Number */
//    sliceThickness = t->sliceThickness;
    dectorYoffset = t->dectorYoffset;
    dectorZoffset = t->dectorZoffset;
    XOffSet = t->XOffSet;
    YOffSet = t->YOffSet;
    ZOffSet = t->ZOffSet;
    phantomXOffSet = t->phantomXOffSet;
    phantomYOffSet = t->phantomYOffSet;
    phantomZOffSet = t->phantomZOffSet;
    XN = RecSize;
    XNC = (XN-1)*0.5;
    YN = RecSize;
    YNC = (YN-1)*0.5;
    ZN = RecSizeZ;
    ZNC = (ZN-1)*0.5;
	YLC    = (YL-1)*0.5+dectorYoffset;
	ZLC    = (ZL-1)*0.5+dectorZoffset;
	dYL= DecWidth/YL;

	dZL= DecHeight/(ZL);       /*Each move Angle of projection*/
	DeltaFai = 2*pi/N_2pi;
    N_pi = N_2pi/2;

    dx = 2*ObjR/XN;
	dy = 2*ObjR/YN;
    dz = 2*ObjR/ZN;

     //////begin of the  main code
    float x,y,z,Dey,Dez,touying,UU,U1,V1,Beta0,Yr,Zr,View,weight;
    int ProjInd,xi,yi,zi,U,V,s0,s1,s2,d1,d2;

	 for(zi = 0; zi<ZN; zi++)
	 {
		 ///compute the projection position for every grid on the image plane
         z = (zi-ZNC) * dz-ZOffSet-phantomZOffSet;
         Beta0 = 2 * pi * z / h;
         s0 = ceil((Beta0 -BetaS) / DeltaFai);
         s1 = s0-N_pi;
         s2 = s0+N_pi-1;

         if ((s1<PN)||(s2>0))
         {
           if (s1 < 0)  {s1 = 0;}
           if (s2 > PN-1) {s2 = PN-1;}
         for (ProjInd = s1; ProjInd <= s2; ProjInd++ )
         {
             View = BetaS + ProjInd* DeltaFai;
             for(yi=0;yi<YN;yi++)
             {
                 y = (yi-YNC)*dy-YOffSet-phantomYOffSet;
                 for(xi=0;xi<XN;xi++)
                 {
                    x  = (xi-XNC)*dx-XOffSet-phantomXOffSet;
                    UU = x*cos(View)+y*sin(View);
				    Yr = x*sin(View)-y*cos(View);
                    Zr = (z-h*View/(2.0*pi))*((DistD)*ScanR)/(UU*sqrt(ScanR*ScanR-Yr*Yr)+ScanR*ScanR-Yr*Yr);
                    U1 = Yr/dYL+YLC;
                    U  = ceil(U1);
                    V1 = Zr/dZL+ZLC;
                    V  = ceil(V1);
                    Dey = U-U1;
                    Dez = V-V1;

                    if ((U>0) && (U<YL) && (V>0) && (V<ZL))
                    {
                           touying = Dey*Dez*t->GF[U-1][V-1][ProjInd]
                                    +Dey*(1-Dez)*t->GF[U-1][V][ProjInd]
                                   +(1-Dey)*Dez*t->GF[U][V-1][ProjInd]
                                   +(1-Dey)*(1-Dez)*t->GF[U][V][ProjInd];
                           weight = 0.5;
                           t->RecIm[xi][yi][zi]=t->RecIm[xi][yi][zi]+weight*touying*DeltaFai;
                    }

                 }// for xi
			 }//for yi
         }//ProjInd
         }//if ((s1<PN)||(s2>0))
	 } //zi
     //////end of the main code
 }

