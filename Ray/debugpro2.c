#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sim5lib2.h"
#include "grishInclude.h"

#define norm_const  2.162000e-45
#define thomson_cross   6.652400e-25 
#define rg2cm       1.485228e+05

static const stokesp_rad EmptyStruct;
    /* const int image_dim_x = grishNx;//2000;
    const int image_dim_y = grishNx;//1500; */
 
    double grishRadialPos; //ALL INCLUDED GLOBALLY HERE AS DEFINED AS EXTERN IN "grishInclude.h" -> (14.09.17) 
    double grishBlobSize;
    int grishB;
    int grishRg;
    double grishVr;
    int grishImg;
    double grishVz;
    double grishStartAngle;
    int grishStartFrame;
    double grishVis;
    double grishExp;
    int grishFrames;
    int gThickdisk;
    int grishModel;
        double grishTime;
    double grishA;            // BH spin (0..1)
        double grishM;           // BH mass (in solar masses)
    double grishDist;            // BH distance (in parsecs)
        double grishIncl;   // observer's inclination (in radians)
        double grishRmax;           // maximal domain radius (in GM/c2 units; max impact parameter to take)
        double grishStep;            // step along the geodesics (in GM/c2 units; only in case of 3D)
        int    grishNx;              //100;            // image resolution (# of pixels over 2*rmax)
//      double param_t    = (argc>1) ? atof(argv[1]) : hhh;  //0.0; // time slice to calculate (0..1 - in fraction of Keplerian                                                 orbital time of the spot)
     
void magspot_main(double a, double M, double dist, double inc, double rmax, double step, double t, int param_nx, stokesp_rad **image, long N_steps, long N_err);


    
    stokesp_rad sp_r; 
    //total_sp.i_r = total_sp.u_r = total_sp.q_r = 0.0;
    //stokesp_rad **image;      //commented out by munawwar 14.02.2018
    
    
//  const double xx         = (pow((v_max/v_min), (1/NN)) - 1) * NN;

 
int main(int argc, char *argv[])
//***************************************************
{   
     

    

    fscanf(stdin,"%d %d %d %d %d %d %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf",&grishImg, &grishNx, &grishFrames, &grishB, &grishModel, &gThickdisk, &grishA, &grishIncl,
    &grishRg, &grishRadialPos, &grishBlobSize, &grishStep, &grishM, &grishDist, &grishRmax, &grishTime, &grishStartAngle, &grishStartFrame, &grishVz, &grishVr, &grishExp, &grishVis);

    double hhh=grishTime;
    
    double param_a    = grishA;            // BH spin (0..1)
        double param_M    = grishM;           // BH mass (in solar masses)
    double param_dist = grishDist;            // BH distance (in parsecs)
        double param_inc  = deg2rad(grishIncl);   // observer's inclination (in radians)
        double param_rmax = grishRmax;           // maximal domain radius (in GM/c2 units; max impact parameter to take)
        double param_step = grishStep;            // step along the geodesics (in GM/c2 units; only in case of 3D)
        int    param_nx   = grishNx;              //100;            // image resolution (# of pixels over 2*rmax)
//      double param_t    = (argc>1) ? atof(argv[1]) : hhh;  //0.0; // time slice to calculate (0..1 - in fraction of Keplerian orbital time of the spot)
    int    param_3D   = grishImg;              // 2D or 3D (0=2D, 1=3D)
    int    param_frames = grishFrames; //80in2d 40in3d         // number of time frames
    double r0           = grishRadialPos; //= grishRadialPos; //*isco(param_a);//((double)(r_ms(param_a)));// <--     //initial radial                                      position of the spot; please introduce by hand !!!
//  sp_r.N_steps=0, sp_r.N_err=0;



if(param_3D == 0) 

{


/*  if (grishRg) {
        r0 = grishRadialPos;
     } else {
        r0 = grishRadialPos*isco(param_a);
     }
*/  


//      int    NN     = 300;            // Number of Bins for energies

//-----------------------------------------Start of time-evolution - for LIGHT CURVE ->

//    double param_tmax = 2.0*M_PI*(sqrt(r0*r0*r0*isco(param_a)*isco(param_a)*isco(param_a))+param_a); // Grish
    double param_tmax = 2.0*M_PI*(sqrt(r0*r0*r0)+param_a); // Grish

//double param_tmax = 2.0*M_PI*(sqrt(r0*r0*r0*6.*6.*6.)+param_a); //1 orbital period //a0->6, a0.5->4.233(4.24)a=0.99->1.4545
    //double param_tmax = 2.0*M_PI*(sqrt(r0*r0*r0*4.2330025*4.2330025*4.2330025)+param_a); //1 orbital period //a0->6, a0.5->4.2330025 a=0.98->  1.6140296
    //double param_tmax = 2.0*M_PI*(sqrt(r0*r0*r0*1.6140296*1.6140296*1.6140296)+param_a); //1 orbital period //a0->6, a0.5->4.233(4.24)a=0.99->1.4545
    //  double param_nu   = 200.0;
   

     
     /// for (a=0.01; a<=0.99; a+=0.49) { 
    ////for (param_inc=deg2rad(22.5); param_inc<=deg2rad(67.5); param_inc+=deg2rad(22.5)) {
    ///param_a = a;
    
    double t=0.0;
    
    // go over time frames
    printf("#a=%.3f incl=%.2f r0=%.3f Tmax=%.3f Spotsize=%.3f\n", param_a, rad2deg(param_inc), r0, param_tmax, grishBlobSize);
    //printf("#a incl r0\n");
    //printf("%.3f %.2f %.3f\n", param_a, param_inc, r0);
    printf("#time\tdeg\tang\ti\tu/i\tq/i\tstep\terr\n"); 
    printf("#---------------\n");
    
    //  printf("#time i ang deg\n"); //<-mine
        fprintf(stderr, "param_frames %.3f\n",(double)(param_frames));



    //for (t=0.0; t<=param_tmax; t+=param_tmax/((double)(param_frames))) {  // original by Michal
    for (t=grishStartFrame*param_tmax/(double)(param_frames); t<=param_tmax; t+=param_tmax/((double)(param_frames))) {    // modified by me, G.

        fprintf(stderr, "frame %e\n", t);


    // allocate image
        stokesp_rad **image;
        image = (stokesp_rad**)calloc(param_nx, sizeof(stokesp_rad*));
        for (int i=0; i<param_nx; i++) image[i] = (stokesp_rad*)calloc(param_nx, sizeof(stokesp_rad));


        // do the calculation 
        magspot_main(param_a, param_M, param_dist, param_inc, param_rmax, param_step, t, param_nx, image, sp_r.N_steps, sp_r.N_err);


//SAVE IMAGE HERE (find how to do) ---> ....

    //for (int j=0; j<param_nx; j++) for (int i=0; i<param_nx; i++)
    //printf("%d %d %.3e %.3e %.3e %.3f\n", i, j, image[j][i].i_r, image[j][i].u_r, image[j][i].q_r, angle);            //Added on 25 Jan 2018


//PRINT IMAGE OR LIGHT-CURVE BELOW>

//------------ sum the result FOR LIGHTCURVE-------------(06.12.17)---------part from main()-----------------------------------------------------
    
        sp_r.i_r = sp_r.u_r = sp_r.q_r = 0.0;
        double maxima_i = image[0][0].i_r;
        for (int i=0; i<param_nx; i++) for (int j=0; j<param_nx; j++) {
            if (maxima_i < image[j][i].i_r) maxima_i = image[j][i].i_r;
        sp_r.i_r += image[j][i].i_r;
            sp_r.u_r += image[j][i].u_r;
            sp_r.q_r += image[j][i].q_r;
        }
        

//SAVE (sp_r)'s HERE (find how to do) ---> ....

        for (int i=0; i<param_nx; i++) free(image[i]);
        free(image);

        // calculate and print the total polarization angle and degree
        double ang = (sp_r.i_r>0.0) ? 0.5*atan2(sp_r.u_r/sp_r.i_r,sp_r.q_r/sp_r.i_r) : 0.0;
        //    double deg = (sp.i>0.0) ? sqrt(sqr(sp.q)+sqr(sp.u))/sp.i : 0.0;
        double deg = (sp_r.i_r>0.0) ? sqrt((sp_r.q_r)*(sp_r.q_r) + (sp_r.u_r)*(sp_r.u_r))/sp_r.i_r : 0.0;
        while (ang < 0.0)  ang += 2.*M_PI;
        while (ang > M_PI) ang -= M_PI;
        printf("%.3f\t%.4f\t%.2f\t%.3e\t%.3e\t%.3e\t%li\t%li\n", t, deg*100., rad2deg(ang), sp_r.i_r, (sp_r.i_r>0.0)?sp_r.u_r/sp_r.i_r:0.0, (sp_r.i_r>0.0)?sp_r.q_r/sp_r.i_r:0.0, sp_r.N_steps, sp_r.N_err); //<-michal
        //printf("%.3f %.3e %.2f %.4f\n", t, sp.i, rad2deg(ang), deg*100.); //<-mine
        //printf("%.3f %.3e\n", t, sp.i);//<-mine
        
        fflush(stdout);
 
}           //TIME STEP FOR LOOP ENDS
//free(&sp_r);

     /* total_sp.i_r += sp_r.i_r;
        total_sp.u_r += sp_r.u_r;
        total_sp.q_r += sp_r.q_r; */


    //if(error) continue;                               //Added by Munawwar 06.02.2018

}           //IF STATEMENT TRUE CONDITION ENDS
else{
    double param_tmax = 2.0*M_PI*(sqrt(r0*r0*r0)+param_a); // Grish

    printf("#a=%.3f incl=%.2f r0=%.3f Tmax=%.3f Spotsize=%.3f\n", param_a, rad2deg(param_inc), r0, param_tmax, grishBlobSize);

    printf("#pix_i\tpix_j\tp_deg\tp_ang\ti\tu/i\tq/i\tfreq\tstep\terr\n");
    printf("#---------------\n");


 // allocate image
    stokesp_rad **image;
    image = (stokesp_rad**)calloc(param_nx, sizeof(stokesp_rad*));
    for (int i=0; i<param_nx; i++) image[i] = (stokesp_rad*)calloc(param_nx, sizeof(stokesp_rad));

    // do the calculation
    magspot_main(param_a, param_M, param_dist, param_inc, param_rmax, param_step, hhh, param_nx, image, sp_r.N_steps, sp_r.N_err);
//fprintf(stderr, "i) DO YOU KNOW\n\tN_step=%li, N_err=%li\n", sp_r.N_steps, sp_r.N_err);

    // print the result
    
     double sp_ang, sp_deg;
     
    printf("#Resolution is %d by %d\n", param_nx, param_nx);

    for (int i=0; i<param_nx; i++) for (int j=0; j<param_nx; j++) {
        stokesp_rad sp_img = image[j][i];
        sp_ang = (sp_img.i_r>0.0) ? 0.5*atan2(sp_img.u_r/sp_img.i_r,sp_img.q_r/sp_img.i_r) : 0.0;
        //  sp_deg = (sp.i>0.0) ? sqrt(sqr(sp.u)+sqr(sp.q))/sp.i : 0.0;
        sp_deg = (sp_img.i_r>0.0) ? sqrt((sp_img.u_r)*(sp_img.u_r) + (sp_img.q_r)*(sp_img.q_r))/sp_img.i_r : 0.0;
        while (sp_ang < 0.0)  sp_ang += 2.*M_PI;
        while (sp_ang > M_PI) sp_ang -= M_PI;
        printf(
            "%d\t%d\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%li\t%li\n", i, j, sp_deg*100., rad2deg(sp_ang), (sp_img.i_r>0)?sp_img.i_r:0.0, (sp_img.i_r>0.)?sp_img.u_r/sp_img.i_r:0.0, (sp_img.i_r>0.)?sp_img.q_r/sp_img.i_r:0.0, sp_img.vu, image[i][j].N_steps, sp_img.N_err);

                                } //for loop ends

for (int i=0; i<param_nx; i++) free(image[i]);
    free(image);
    printf("\n");

}                       //IF() for IMAGE ends here


//printf("%d %d %.3e %.3e %.3e %.3e %.3li %.3li\n", i, j, image[j][i].i_r, image[j][i].u_r/image[j][i].i_r, image[j][i].q_r/image[j][i].i_r, image[j][i].p_ang, N_err, N_steps); //COMMENTED OUT by Munawwar on 14.02.2018



    return 0;
}



void magspot_main(double a, double M, double dist, double inc, double rmax, double step, double t, int param_nx, stokesp_rad **image, long N_steps, long N_err)
// main routine for calculations
// a      ... BH spin (0..1)
// M      ... BH mass (in solar mass units)
// dist   ... BH distance (in parsecs)
// i      ... observer's inclination (0..85)
// rmax   ... maximal domain radius (in GM/c2 units; max impact parameter to take)
// step   ... step along the geodesics (in GM/c2 units; only in case of 3D)gRadi
// time   ... time slice to calculate (0..1 - in fraction of Keplerian orbital time of the spot)
// calc3D ... 2D or 3D (0=2D, 1=3D)
// N      ... size of the image
// image  ... image (the result)


{


//  void image_alpha_beta(double alpha, double beta, double dOmega,double a, double M, double inc, double rmax, double step, double t, int calc3D, stokesp_rad *sp_r);

//                          Ray-tracing code beginning


 /*   // read paramters from command line
    double a   = grishA; //(argc==3) ? atof(argv[1]) : 0.0;
    double inc = grishIncl; //(argc==3) ? deg2rad(atof(argv[2])) : deg2rad(65.); */
    
  /* if ((a<0.0) || (a>0.999) || (inc<0.0) || (inc>deg2rad(89.))) {
        fprintf(stderr, "ERROR: parameters out of range\n");
        exit(0);
    }
*/

 /*   // allocate memory for images (2D arrays of float)
    // we have arrays for flux, g-factor
    float* image = (float*)calloc(image_dim_x*image_dim_y, sizeof(float));   */

    // width of view - how big portion of the disk will the image cover
    double rms  = r_ms(a);
    //double rbh  = r_bh(a);
    double view_width = rms + 8.0;
    
    //float precision = 1.0;
    //float r_sphere = (float)(grishBlobSize);
    float r_max = (float)(rmax);
        
    // start the clock and calculate image
 /*   clock_t cpu_t1, cpu_t2;
    cpu_t1 = clock(); */



/*    int ix, iy;
    // go over image and for each pixel independently determine its brightness
    for (iy=0; iy<image_dim_y; iy++) for (ix=0; ix<image_dim_x; ix++) { */
 //       if (ix==0) fprintf(stderr,"y=%d\n",iy);
        // impact parameters
        // these are image coordinates scaled to view_width and shifted so 
        // that [0,0] is in the middle of the image

    image[0][2].i_r = 45.000; 
    image[0][0].i_r = 3333.000; 
    
    for (int j=0; j<param_nx; j++) for (int i=0; i<param_nx; i++) {

        double alpha = (((double)(i)+.5)/(double)(param_nx)-0.5)*2.0*view_width;
        double beta  = (((double)(j)+.5)/(double)(param_nx)-0.5)*2.0*view_width * ((double)param_nx/(double)param_nx);
    
        N_steps=0; N_err=0;
    //stokesp_rad *sp_r; sp_r = (struct stokesp_rad *) malloc (sizeof(struct stokesp_rad));             //Added by Munawwar 01.02.2018
        int error;
        double f;double x[4];
        geodesic g, gd, debuga, trace;






 // initialize photon trajectory for particular impact parameters [alpha,beta]
        
     int check = 1; 

    gd = geodesic_init_inf_new(inc, a, alpha, beta, &g, &error, check, INFINITY);


//    fprintf(stderr, "a=%e, inc = %e, check value is %d and error is %d, g.Rpc=%e, gd.Rpc=%e\n", a, inc, check, &error,g.Rpc,gd.Rpc);
     if (check != 1) {
        fprintf(stderr, "photon rejected with error %d %e\n", error, gd.q);
        continue;  
    }
 
    
    //fprintf(stderr, "To see for TRACE: Rpc=%e, rp=%e\n", trace->Rpc, trace->rp);
    double Rpc=gd.Rpc, rp=gd.rp;
    double P = geodesic_P_int_new(gd, r_max, 0);
    if(Rpc < rp) { 
//        fprintf(stderr, "(I) Check values!\n\tr_max = %f , P = %e , 2Rpc = %e and debug-rp=%e, trace_rp=%e, cr1=%e+%ei, cr2=%e+%ei, cr3=%e+%ei, cr4=%e+%ei\n",
//            r_max, P, 2.*Rpc, debuga.rp, trace.rp, creal(gd.r1), cimag(gd.r1), creal(gd.r2), cimag(gd.r2), creal(gd.r3), cimag(gd.r3), creal(gd.r4), cimag(gd.r4));
                }

 /*        double r1 = creal(gd.r1);
         double r2 = creal(gd.r2);
         double r3 = creal(gd.r3);
         double r4 = creal(gd.r4);
         double mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
         double RPP  = r1;
         double RPCC = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt((r2-r4)/(r1-r4)), mm);

         double m4 = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
         double x4 = 0.5*fabs(P - RPCC)*sqrt((r2-r4)*(r1-r3));
         double sn2 = pow(jacobi_sn(x4,m4), 2.0);
         double R_main = ( r1*(r2-r4)-r2*(r1-r4)*sn2 ) / ( r2-r4-(r1-r4)*sn2 );

fprintf(stderr, "(II) CHECK THIS -\n\tr1=%e,r2=%e,r3=%e,r4=%e,mm=%e, RPP=%e, RPCC=%e, m4=%e, x4=%e, sn2=%e, MAINVALUE-R=%e\n", r1,r2,r3,r4,mm,RPP,RPCC,m4,x4,sn2,R_main);
*/


    x[0] = 0.0;  // zero relative time delay
    x[1] = geodesic_position_rad_new(gd, P);   // this is unnecessary as it is supposed to get r=1000.0
    x[2] = geodesic_position_pol_new(gd, P);
    x[3] = 0.0;  // zero initial phi angle
    
//fprintf(stderr, "(III) Check again!\n\tx[1] = %e , x[2] = %e,  r_max = %f , P = %e , 2Rpc = %e and rp = %e\n", x[1], x[2], r_max, P, 2.*Rpc, rp);

    stokesp_rad model(geodesic gd, double t, double r, double m, double phi, double inc, double step, int *error, double alpha, double beta, long N_err, long N_steps, double P);

    image[j][i] = model(gd, x[0], x[1], x[2], x[3], inc, step, &error, alpha, beta, N_err, N_steps, P);
//  fprintf(stderr, "ii) DO YOU KNOW\n\tN_step=%li, N_err=%li, intensity=%e\n", image[j][i].N_steps, image[j][i].N_err, image[j][i].i_r);

    
    /*  image[j][i].i_r=sp_r.i_r;
    image[j][i].u_r=sp_r.u_r;
    image[j][i].q_r=sp_r.q_r;
    */

       
    } //for loop of image calculation ends here 

//fprintf(stderr, "Intensity at [0,2]:%e\tSteps in [0,2]:%li\tFreq=%e\n", image[0][2].i_r, image[0][2].N_steps, image[0][2].vu);
//image[0][0].i_r=3.04595; image[0][0].N_steps=54;

return;

//--------------------------------------------------------------------------------------------------------------------------------------------------
    // end if image calculation, read the clock
/*    cpu_t2 = clock();

    fprintf(stderr,
        "profiling: N=%d t=%.2es rate=%.3e per s\n", 
        param_nx*param_nx, (double)(cpu_t2-cpu_t1)/CLOCKS_PER_SEC, 
        param_nx*param_nx/((double)(cpu_t2-cpu_t1)/CLOCKS_PER_SEC)
    ); */

//--------------------------------------------------------------------------------------------------------------------------------------------------

}

//-----------------------------------------------------Ray-tracing code end-------------------------------------------------------




stokesp_rad model(geodesic gd, double t, double r, double m, double phi, double inc, double step, int *error, double alpha, double beta, long N_err, long N_steps, double P)

    
    {

    int RadTrans = 1;
    double a=gd.a;
    // find where the trajectory crosses radius r=1000.
    // obtain position parameter (returns NaN if there is no crossing)
    float r_max = (float)(grishRmax);
    double rbh  = r_bh(a);
    float precision = 1.0;
    float r_sphere = (float)(grishBlobSize);
    double mu, vu_loc, Kv, Tau=0.0, del_Tau;
    double x[4];
    double g_factor;
    double Rho, B_const, D_i;
    double dl;

    stokesp_rad total_sp;
    total_sp.N_steps=N_steps; 
    total_sp.N_err=N_err;
    double min_T = 1e-99; 
    total_sp.vu=2.3e11; 
    total_sp.T_factor=1.0;
    total_sp.i_r=0.0; 
    total_sp.u_r=0.0; 
    total_sp.q_r=0.0;

    x[0] = 0.0;  // zero relative time delay
    x[1] = geodesic_position_rad_new(gd, P);   // this is unnecessary as it is supposed to get r=1000.0
    x[2] = geodesic_position_pol_new(gd, P);
    x[3] = 0.0;  // zero initial phi angle
    

    double k[4];
    geodesic_momentum_new(gd, P, x[1], x[2], k);

    /*double dm = geodesic_dm_sign_new(gd, P);
    double r_sign=P<gd.Rpc?-1:+1;
    double m_sign=dm;
    double a2 = sqr(gd.a);
    double l2 = sqr(gd.l);
    double r2 = sqr(r);
    double m2 = sqr(m);
    double S = r2 + a2*m2;
    double D = r2 - 2.*r + a2;
    double R = sqr(r2+a2-gd.a*gd.l) - D*(sqr(gd.l-gd.a) + gd.q);
    double M = gd.q - l2*m2/(1.-m2) + a2*m2;
    k[0] = +1/S * ( -gd.a*(gd.a*(1.-m2)-gd.l) + (r2+a2)/D*(r2+a2-gd.a*gd.l) );
    k[1] = +1/S * sqrt(R);
    k[2] = +1/S * sqrt(M);
    k[3] = +1/S * ( -gd.a + gd.l/(1.-m2) + gd.a/D*(r2+a2-gd.a*gd.l));
    if (r_sign<0.0) k[1] = -k[1];
    if (m_sign<0.0) k[2] = -k[2];
 
    fprintf(stderr, "a2=%e,l2=%e,r2=%e,m2=%e,S=%e,D=%e,R=%e,M=%e,k1=%e,k2=%e,k3=%e,k4=%e\n", a2,l2,r2,m2,S,D,R,M,k[0],k[1],k[2],k[3]); */
//    fprintf(stderr, "(IV.I) AGAIN !!\n\tx1=%e, x2=%e, r = %e , r_max = %f , P = %e , 2Rpc = %e and rp = %e; k0=%e, k1=%e, k2=%e, k3=%e\n", x[1], x[2], r, r_max, P, 2*gd.Rpc, gd.rp, k[0],k[1],k[2],k[3]);
    if (isnan(k[0])) 
    {
        fprintf(stderr, "(I) photon rejected with k=NAN (%e)\n", gd.q); 
        return EmptyStruct;                                             //RETURN EmptyStruct added 06.02.2018
        //continue;
    }

    
    //total_sp = EmptyStruct;
    //total_sp = (struct stokesp_rad *) malloc (sizeof(struct stokesp_rad));                //Added by Munawwar 01.02.2018

    raytrace_data rtd;
    raytrace_prepare(a, x, k, precision, 0, &rtd);
 //   fprintf(stderr, "(IV.II) Once for k[0] ->\n\tk0=%e, k1=%e, k2=%e, k3=%e\n\tx0=%e, x1=%e, x2=%e, x3=%e\n",k[0],k[1],k[2],k[3],x[0],x[1],x[2],x[3]);

    dl=1.0; //First time while is called, dl=1.0 and then changes accordingly during while loop.

    while(1)
    {
    //    dl = (x[1] < r_sphere) ? 5e-2  : 1.0; // limited max step in the sphere
    //    fprintf(stderr, "dl=%e\n", dl);

        raytrace(x, k, &dl, &rtd);
    //    fprintf(stderr, "(IV.III) For k[0] ->\n\tk0=%e, k1=%e, k2=%e, k3=%e\n\tx0=%e, x1=%e, x2=%e, x3=%e\n",k[0],k[1],k[2],k[3],x[0],x[1],x[2],x[3]);
        if (isnan(k[0])) 
        {

            fprintf(stderr, "(II) photon rejected with k=NAN (%e)\n", gd.q); 
            //total_sp=EmptyStruct;
            break;                                            //RETURN EmptyStruct added 06.02.2018
            //continue;
        }

        // if (dl < 1e-3) fprintf(stderr,"dl=%e,  a=%e  b=%e  r=%e\n", dl, alpha, beta, x[1]);
               
        // stop condition:
        // photon either reaches black hole or escapes beyond r_max bounday
        if ((x[1] < rbh*1.05) || (x[1] > r_max*1.1))
        {
            //fprintf(stderr, "Black Hole reached or photon escaped beyond R_max\n"); 
            //total_sp=EmptyStruct;
            total_sp.N_steps++;
            break;
        }

        // stop if relative error in this step is too large
        if (rtd.error>1e-1) 
        {
            fprintf(stderr, "raytrace(): aborted due to large error at r=%e m=%e\n", x[1], x[2]);
            total_sp.N_err += 1;
            total_sp.N_steps ++;
            continue;
            //total_sp=EmptyStruct;
            //break;
        }

        // increase accumulated emissivity while inside the sphere
        // Iv += j * ds
        /*if (r < r_sphere) 
        {*/
    //            fprintf(stderr, "(IV.IV.1) For k[0] ->\n\tk0=%e, k1=%e, k2=%e, k3=%e\n\tx0=%e, x1=%e, x2=%e, x3=%e\n",k[0],k[1],k[2],k[3],x[0],x[1],x[2],x[3]);
       
       sim5metric metric;
     //  metric = (sim5metric *) malloc(sizeof(sim5metric));
       
 
       
       t = x[0]; r=x[1]; m = x[2]; phi = x[3];      
       double k_0 = k[0], k_1 = k[1], k_2 = k[2], k_3 = k[3];




       kerr_metric(a, r, m, &metric);
        
       /*
        double r22  = x[1]*x[1];
        double a22  = sqr(a);

        double m22  = x[2]*x[2];
        double SS   = r22 + a22*m22;
        double s2S = (1.0-m22)/SS;
        metric->a = a;
        metric->r = x[1];
        metric->m = x[2];
        metric->g00 = -1. + 2.0*(x[1]/SS);
        metric->g11 = SS/(r22-(2.*x[1])+a22); 
        metric->g22 = SS;
        metric->g33 = (r22*r22 + a22*a22*m22 + a22*r22*(1.+m22) + 2.*x[1]*a22*s2S*SS)*s2S;
        metric->g03 = -2.*a*x[1]*s2S;
        */




    //   fprintf(stderr, "(IV.IV.2) For k[0] ->\n\tk0=%e, k1=%e, k2=%e, k3=%e\n\tx0=%e, x1=%e, x2=%e, x3=%e\n",k[0],k[1],k[2],k[3],x[0],x[1],x[2],x[3]);
       double Omega = OmegaK(r,a);                         // Omega of uniform Keplerian rotation --> no shear

        sim5tetrad tetrad;
     //   tetrad = (sim5tetrad*) malloc(sizeof(sim5tetrad));

         tetrad_azimuthal(&metric, Omega, &tetrad);

     /*         if (Omega==0.0) 
       {
         tetrad->e[0][0] = sqrt(metric->g33/(sqr(metric->g03) - metric->g33*metric->g00));
         tetrad->e[0][1] = 0.0;
         tetrad->e[0][2] = 0.0;
         tetrad->e[0][3] = -tetrad->e[0][0] * metric->g03/metric->g33;

         tetrad->e[1][0] = 0.0;
         tetrad->e[1][1] = 1./sqrt(metric->g11);
         tetrad->e[1][2] = 0.0;
         tetrad->e[1][3] = 0.0;

         tetrad->e[2][0] = 0.0;
         tetrad->e[2][1] = 0.0;
         tetrad->e[2][2] = -1./sqrt(metric->g22);
         tetrad->e[2][3] = 0.0;

         tetrad->e[3][0] = 0.0;
         tetrad->e[3][1] = 0.0;
         tetrad->e[3][2] = 0.0;
         tetrad->e[3][3] = 1./sqrt(metric->g33);
         tetrad->metric = *metric;
       }

                double g00 = metric->g00;
                double g33 = metric->g33;
                double g03 = metric->g03;
                double U0 = sqrt(-1.0/(g00 + 2.*Omega*g03 + sqr(Omega)*g33));
                double U3 = U0*Omega;

                tetrad->e[0][0] = U0;
                tetrad->e[0][1] = 0.0;
                tetrad->e[0][2] = 0.0;
                tetrad->e[0][3] = U3;

                tetrad->e[1][0] = 0.0;
                tetrad->e[1][1] = sqrt(1./metric->g11);
                tetrad->e[1][2] = 0.0;
                tetrad->e[1][3] = 0.0;

                tetrad->e[2][0] = 0.0;
                tetrad->e[2][1] = 0.0;
                tetrad->e[2][2] = -sqrt(1./metric->g22);  
                tetrad->e[2][3] = 0.0;

                double k1 = (g03*U3+g00*U0);
                double k2 = (g33*U3+g03*U0);
                tetrad->e[3][0] =  - sign(k1)*k2 / sqrt((g33*g00-g03*g03)*(g00*U0*U0+g33*U3*U3+2.0*g03*U0*U3));
                tetrad->e[3][1] = 0.0;
                tetrad->e[3][2] = 0.0;
                tetrad->e[3][3] = tetrad->e[3][0] * (-k1/k2);

                tetrad->metric = *(metric);
        */

       double U[4];
       //Next 3 lines taken from sim5lib.c to declare the vaf0=x-x0-1.5*a*log(x/x0);riables (08.09.17)  
       double gtt = -1. + 2./r;
       double gtf = -2. * a/r;
       double gff = sqr(r) + sqr(a) + 2. * sqr(a)/r;
       U[0] = sqrt(-1.0/(gtt + 2.*Omega*gtf + sqr(Omega)*gff));
       U[1] = 0.0;
       U[2] = 0.0;
       U[3] = U[0]*Omega;
       double disk_nt_disk_rms = 6.0;
       double xx=sqrt(r);
       double x0=sqrt(disk_nt_disk_rms);
       double f0=xx-x0-1.5*a*log(xx/x0);

       double phi0 = f0 + Omega*x0; //angular position of spot - changed x[0] to x0 - 22.01.2018
       //   double phi0 = params.F0 + Omega0*t ............OLDER STATEMENT GIVES ERRORS

       double dX[4];
       dX[0]=0; 
       dX[1] = r-grishRadialPos; 
       dX[2] = acos(m)-(M_PI/2); 
       dX[3] = phi - phi0;

       double p_t = k[0]*gtt + k[3]*gtf;
    //   fprintf(stderr, "FURTHER:\n\tgtt=%e, gtf=%e, gff=%e\n\t, k[0]=%e, k[1]=%e, k[2]=%e, k[3]=%e\n\t x[0]=%e, x[1]=%e, x[2]=%e, x[3]=%e\n", gtt, gtf, gff, k[0], k[1], k[2], k[3], x[0], x[1], x[2], x[3]);

       double p_x_u = k[0]*U[0]*gtt + metric.g11*k[1]*U[1] + metric.g22*k[2]*U[2] + k[3]*U[3]*gff + k[0]*U[3]*gtf + k[3]*U[0]*gtf; //08.09.17 ---- changed from grr to g11 and ghh to g22
       //double p_x_u = p_t/2.;

       g_factor = p_t/p_x_u;    

    //   fprintf(stderr, "G_factor=%.3e\n", g_factor);

    //   fprintf(stderr, "(IV.IV.3) For k[0] ->\n\tk0=%e, k1=%e, k2=%e, k3=%e\n\tt=%e, r=%e, m=%e, phi=%e\n\tgtt=%e, gtf=%e, gff=%e, x=%e\n\tx0=%e, f0=%e, phi0=%e, Omega=%e, U0=%.3e\n",k[0],k[1],k[2],k[3],t,r,m,phi,gtt,gtf,gff,xx,x0,f0,phi0,Omega,U[0]);
    
       double R2 = dotprod(dX,dX, &metric) + pow(dotprod(dX,U, &metric),2);

       dl = (sqrt(R2) < r_sphere) ? 5e-2  : 1.0;

       double emis;
       double K_em = exp(-0.5*R2/sqr(grishBlobSize));
       if(K_em < exp(-grishVis)) emis=0;
       emis=K_em*0.03/pow(grishBlobSize,3)*dl;


       int indicator = 1;
       indicator = (R2 > pow(grishBlobSize,2.) || isnan(R2)) ? 1 : 0 ;
    //   fprintf(stderr, "(V) (Indicator)\n\tIndicator=%d, R(squared)=%e, Blob(squared)=%e, r=%e, RadialPos=%e\n", indicator, R2, pow(grishBlobSize,2),r,grishRadialPos);

    //   fprintf(stderr, "(IV.IV.4) For k[0] ->\n\tk0=%e, k1=%e, k2=%e, k3=%e\n\tx0=%e, x1=%e, x2=%e, x3=%e\n\tgtt=%e, gtf=%e, gff=%e, x=%e\n\tx-0=%e, f0=%e, phi0=%e, Omega=%e\n",k[0],k[1],k[2],k[3],x[0],x[1],x[2],x[3],gtt,gtf,gff,xx,x0,f0,phi0,Omega);

       //  fprintf(stderr, "i This actually runs! MY GOODNESS\n");
      

     




   //fprintf(stderr, "Values Metric:\n\tG00=%e\tG11=%e\tG22=%e\n\tG33=%e\tG03=%e\tU0=%e\n", metric.g00,metric.g11,metric.g22,metric.g33,metric.g03, sqrt(-1.0/(metric.g00 + 2.*Omega*metric.g03 + sqr(Omega)*metric.g33))); 
    
    //   fprintf(stderr, "(IV.IV.5) For OLD k[0] ->\n\tk0=%e, k1=%e, k2=%e, k3=%e\n",k0,k1,k2,k3);

       double B[4];
       //   Declaration of P[] done in the following from sim5lib.h

       //k[0]=k_0;k[1]=k_1;k[2]=k_2;k[3]=k_3;
       //x[0]=t;x[1]=r;x[2]=m;x[3]=phi;
       
    //   fprintf(stderr, "(IV.IV.6) For k[0] ->\n\tk0=%e, k1=%e, k2=%e, k3=%e\n\tx0=%e, x1=%e, x2=%e, x3=%e\n",k[0],k[1],k[2],k[3],x[0],x[1],x[2],x[3]);

       double B_loc[4];
       B_loc[0] = 0.0;  // B.U=0 (B is space-like) => B_loc[0]=0 
       B_loc[1] = 0.0;
       if ( grishB == 0 ) 
        { 
            B_loc[2] = 1.0;  // locally vertical B-field 
            B_loc[3] = 0.0; 
        } 
        else 
        {   
            B_loc[2] = 0.0;  // d/d(theta) aligned field
            B_loc[3] = 1.0;    //locally azimuthal B-field
        }
       // fprintf(stderr, "iii This actually runs! MY GOODNESS\n");

       on2bl(B_loc, B, &tetrad);
        /*int i,j;
        for (i=0;i<4;i++) {
        B[i] = 0.0;
        for (j=0;j<4;j++) B[i] += B_loc[j] * (tetrad->e[j][i]);} */
    //  fprintf(stderr, "B_vect=>\n\tB0=%e\tB1=%e\tB2=%e\tB3=%e\n", B[0],B[1],B[2],B[3]);

       double pol_deg;
       stokesp_rad emission_stokes(double P[4], double U[4], double B[4], sim5metric metric, sim5tetrad tetrad, double alpha, double beta, double emis, stokesp_rad total_sp, double *pol_deg);
       stokesp_rad sp_loc;
        
       sp_loc = emission_stokes(k, U, B, metric, tetrad, alpha, beta, emis, total_sp, &pol_deg);
       //   fprintf(stderr, "%.3f\n" , angle);
      
       //    stokesp_rad sp_r = emission_stokes(g, r, m, F, U, B, poldeg);        //,ds);

       if(indicator==0)            //Added by Munawwar 02.03.18
        {   
        if(RadTrans==1)
        {
            Rho                 =  1e-7 * (1.0 - R2/sqr(grishBlobSize));
            step                =  dl;
            mu                  =  1.0e+11; //Assuming BLR region for accretion for electrons with 1e+06 cm-3
            vu_loc              =  (total_sp.vu)/g_factor;
            Kv                  =  1./(rg2cm*r_max*Rho) ; //Assumption for Alpha=1*r_max
            // double Em_v      =  (4.*M_PI*2.*10e-45*pow(1000*vu_loc*planck_h,-2.5))/(mu*mass_proton);
            B_const             =  (mu*mass_proton)/(Rho*thomson_cross*r_max*rg2cm);
            del_Tau             =  step/(r_max);
            Tau                 += del_Tau;
            D_i                 =  (norm_const*pow((joule2erg*planck_h*vu_loc*erg2kev),-2.5)*total_sp.T_factor*rg2cm*step*exp(-(del_Tau)/2))/(B_const*thomson_cross*pow(vu_loc,2)); //Added by Munawwar 31.01.2018
            //D_i               =  1./(4.*M_PI)*sp_r->T_factor*Rho*Em_v*exp(-(del_Tau)/2.)*step/(pow(vu_loc,2));
            total_sp.i_r        += D_i * pow(g_factor,4);
            total_sp.q_r        += D_i * cos(2 * sp_loc.p_ang);
            total_sp.u_r        += D_i * sin(2 * sp_loc.p_ang);

            total_sp.T_factor   *= exp(-(del_Tau));
        }
        else
        {
            Rho                 =  1e-7 * (1.0 - R2/sqr(grishBlobSize));
            step                =  dl;
            mu                  =  1.0e+11; //Assuming BLR region for accretion for electrons with 1e+06 cm-3
            vu_loc              =  (total_sp.vu)/g_factor;
            Kv                  =  1./(rg2cm*r_max*Rho) ; //Assumption for Alpha=1*r_max
            // double Em_v      =  (4.*M_PI*2.*10e-45*pow(1000*vu_loc*planck_h,-2.5))/(mu*mass_proton);
            B_const             =  (mu*mass_proton)/(Rho*thomson_cross*r_max*rg2cm);
            del_Tau             =  0;
            Tau                 += del_Tau;
            D_i                 =  (norm_const*pow((joule2erg*planck_h*vu_loc*erg2kev),-2.5)*total_sp.T_factor*rg2cm*step*exp(-(del_Tau)/2))/(B_const*thomson_cross*pow(vu_loc,2)); //Added by Munawwar 31.01.2018
            //D_i               =  1./(4.*M_PI)*sp_r->T_factor*Rho*Em_v*exp(-(del_Tau)/2.)*step/(pow(vu_loc,2));
            total_sp.i_r        += sp_loc.i_r * pow(g_factor,4);
            total_sp.q_r        += sp_loc.i_r * pol_deg * cos(2 * sp_loc.p_ang);
            total_sp.u_r        += sp_loc.i_r * pol_deg * sin(2 * sp_loc.p_ang);
    //        fprintf(stderr, "int_loc=%e\tpoldeg=%e\tPol_ang=%e\n",sp_loc.i_r, pol_deg, sp_loc.p_ang);

        }
            if(min_T > total_sp.T_factor) min_T = total_sp.T_factor; 
        
    
    //        fprintf(stderr, "RadTransEndsNow\n\tp_t=%e , p_x_u=%e , g_factor=%e, Rho=%e, dl=%e, step=%e, r=%e\n", p_t, p_x_u, g_factor, Rho, dl, step,r);
            if ( Tau > 1 || min_T < 1e-99 || m <= 0 || isnan(R2)) 
                {   //fprintf(stderr, "Rad Transfer + Geodesic End ->\n\tOptDepth=%.3e\tAbsCoef=%.3e\tTheta=%.3e\tDegInt=%.3e\n\tMetG11=%.3e\tMetG22=%.3e\tP_t=%.3e\tP_x_u=%.3e\tOmega=%.3e\n",Tau,min_T,m,D_i,metric.g11, metric.g22, p_t, p_x_u, Omega);
                    //fprintf(stderr, "\tdTau=%e\tStep=%e\tRvect=%e\tIntensity=%e\n", del_Tau, step, sqrt(R2), total_sp.i_r);
                    //fprintf(stderr, "\tNSteps=%li\tFreq=%e\tu=%e\tq=%e\tDensity=%e\tAng=%e\n", total_sp.N_steps, total_sp.vu, total_sp.u_r, total_sp.q_r, Rho, total_sp.p_ang);
                    total_sp.N_steps++;
                    break;
                }
            
 

            /*            float g = -metric.g00;
                image[ix+image_dim_x*iy] += pow(g,3)*dl; */

        }
          else

        {
            //fprintf(stderr, "Photon on its way to Blob ->\n\tR_sqr=%.3e\t>\tBlob_sqr=%.3e\n", R2, pow(grishBlobSize,2.));
            total_sp.i_r            +=  0.0; 
            total_sp.q_r            +=  0.0;
            total_sp.u_r            +=  0.0;
            total_sp.T_factor       +=  0.0;
            total_sp.p_ang          +=  0.0; 
            total_sp.N_err          +=  0.0;
        }

     
        total_sp.N_steps++;

    } //WHILE ENDS    

//    fprintf(stderr, "WHILE HAS NOW ENDED and N_steps=%li\n", total_sp.N_steps);
//    total_sp.i_r=489.432; total_sp.q_r=97.89932; total_sp.u_r=-234.782107; 
    //free(tetrad);
    //free(metric);
    return total_sp;
        }





stokesp_rad emission_stokes(double P[4], double U[4], double B[4], sim5metric metric, sim5tetrad tetrad, double alpha, double beta, double emis, stokesp_rad total_sp, double *pol_deg)
//***************************************************
{

    //fprintf(stderr, "v This actually runs! MY GOODNESS\n");
     // local 4-velocity
    double U_loc[4];
    bl2on(U, U_loc, &tetrad);
    //fprintf(stderr, "U_vect=>\n\tU0=%e\tU1=%e\tU2=%e\tU3=%e\n", U_loc[0],U_loc[1],U_loc[2],U_loc[3]);
    // local photon momentum (tetrad components of P[4])
    double P_loc[4];
    bl2on(P, P_loc, &tetrad);
    //fprintf(stderr, "P_vect=>\n\tP0=%e\tP1=%e\tP2=%e\tP3=%e\n", P_loc[0],P_loc[1],P_loc[2],P_loc[3]);
    // arbitrary local random vector (tetrad components)
    double X_loc[4]; 
    double rand_limit = 432.0;
    X_loc[0] = 0.0;
    X_loc[1] = ((float)sim5rand()/(float)(RAND_MAX)) * rand_limit;
    X_loc[2] = ((float)sim5rand()/(float)(RAND_MAX)) * rand_limit;
    X_loc[3] = ((float)sim5rand()/(float)(RAND_MAX)) * rand_limit;

    // arbitrary local vector perpendicular to P_loc (tetrad components)
    double V_loc[4];
    V_loc[0] = 0.0;
    V_loc[1] = P_loc[2]*X_loc[3] - P_loc[3]*X_loc[2];
    V_loc[2] = P_loc[3]*X_loc[1] - P_loc[1]*X_loc[3];
    V_loc[3] = P_loc[1]*X_loc[2] - P_loc[2]*X_loc[1];
    
    // arbitrary local vector perpendicular to P_loc & V_loc (tetrad components)
    double W_loc[4];
    W_loc[0] = 0.0;
    W_loc[1] = P_loc[2]*V_loc[3] - P_loc[3]*V_loc[2];
    W_loc[2] = P_loc[3]*V_loc[1] - P_loc[1]*V_loc[3];
    W_loc[3] = P_loc[1]*V_loc[2] - P_loc[2]*V_loc[1];
    
    // transformation of V_loc to B-L coordinates + normalization V*V=1
    double V[4];
    on2bl(V_loc, V, &tetrad);
    vector_norm_to(V, 1.0, &metric);

    // transformation of W_loc to B-L coordinates + normalization W*W=1
    double W[4];
    on2bl(W_loc, W, &tetrad);
    vector_norm_to(W, 1.0, &metric);
    
    // scalar products
    double UP = dotprod(U, P, &metric);
    double UV = dotprod(U, V, &metric);
    double UW = dotprod(U, W, &metric);
    double BP = dotprod(B, P, &metric);
    double BV = dotprod(B, V, &metric);
    double BW = dotprod(B, W, &metric);
    double VV = dotprod(V, V, &metric);
    double WW = dotprod(W, W, &metric);
    double BB = dotprod(B, B, &metric);
//    fprintf(stderr, "Scalar Products:\n\tUP=%e\tUV=%e\tUW=%e\tBP=%e\n\tBV=%e\tBW=%e\tVV=%e\tWW=%e\tBB=%e\n", UP,UV,UW,BP,BV,BW,VV,WW,BB);
    double nX = (UP*BV - UV*BP)/sqrt(sqr(UP)*BB-sqr(BP))/sqrt(VV);
    double nY = (UP*BW - UW*BP)/sqrt(sqr(UP)*BB-sqr(BP))/sqrt(WW);
//    fprintf(stderr, "N_x=%e\tN_y=%e\n\tTerm1=%e\tTerm2=%e\tTerm3=%e\tTerm4=%e\n", nX, nY, (sqr(UP)*BB-sqr(BP)), sqrt(sqr(UP)*BB-sqr(BP)), VV, sqrt(VV) );
    //checks:
    //if (isnan(nY)&&(emi>0)) fprintf(stderr,"emi=%e b0=%e b1=%e b2=%e b3=%e\n", emi, B[0], B[1], B[2], B[3]);
    //if (isnan(nY)&&(emi>0)) fprintf(stderr,"WW=%e BB=%e BP=%e  UP=%e  X=%e\n", WW, BB, BP, UP, sqr(UP)*BB-sqr(BP));
    double i = emis;
    double q = i * (sqr(nY) - sqr(nX));
    double u = i * nX*nY;
    

//    fprintf(stderr, "Check values q and u from Emission Stokes:\n\tQ=%e\tU=%e\n", q, u);
    //checks:
    //if ((isnan(sp_loc.q*sp_loc.u)) fprintf(stderr, "%e UP=%e BB=%e BP=%e BU=%e UU=%e\n", sqr(UP)*BB-sqr(BP), UP, BB, BP, BU, UU);

    sim5complex kappa = photon_wp_const(P,W, &metric);
    double pol_a_inf = polarization_angle_infty(metric.a, grishIncl, alpha, beta, kappa); //BEFORE WAS metric->a

     double angF, ang, deg;
  
    // polarization degree & angle (angle in range [0..PI])
//  ang  = 0.5*atan2(sp_loc.u,sp_loc.q) + polarization_angle_infty(metric->a, /*inclination of observer*/, alpha, beta, kappa); OLD STATEMENT CHANGED TO BELOW on 08.09.17
    ang  = 0.5*atan2(u,q) + pol_a_inf;
    *pol_deg  = sqrt(sqr(u)+sqr(q))/i; 
//  fprintf(stderr, "int_loc=%e\tpoldeg=%e\n",i, pol_deg);

    total_sp.i_r = i;
    total_sp.q_r = q;
    total_sp.u_r = u;
    total_sp.p_ang = ang;

    return total_sp;
}






