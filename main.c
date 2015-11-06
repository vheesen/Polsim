#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"

long i, j, k, gs;

int disk, halo, disk_odd, halo_odd, halo_n, equipartition, gaussian_bfield, clockwise_north, dynamo, idyn, cone;


double x[207], y[207], z[207], y_ss[207][207][207], z_ss[207][207][207], 
    x_ss[207][207], Q_0[207][207], Q_1[207][207], U_0[207][207], U_1[207][207], Q_2[207][207], 
    U_2[207][207], Q_3[207][207], 
    U_3[207][207], RM[207][207], RM_2[207][207], 
    psi_d[207][207][207], psi_h[207][207][207], delta_phi_1[207][207][207], 
    delta_phi_2[207][207][207], delta_phi_3[207][207][207], map[207][207], PI_1[207][207], PI_2[207][207],
    E_x_d[207][207][207], E_y_d[207][207][207],
    E_z_d[207][207][207],
    E_x_h[207][207][207], E_y_h[207][207][207], E_z_h[207][207][207],
    Et[207][207], delta_phi_o1[207][207], PI_0[207][207];

double delta_x, delta_y, delta_z, E, E_phi, E_R, E_phi_h, E_phi_d, E_R_h, E_R_d,
    E_xs_d, E_ys_d, E_zs_d, E_xs_h, E_ys_h, E_zs_h,
    theta, phi, x_ref, y_ref, z_ref, ne, ne_0, h_e, 
    lambda_1, lambda_2, lambda_3, psi_1_d, psi_2_d, psi_3_d, psi_1_h, psi_2_h, psi_3_h, psi_obs_1, psi_obs_2, psi_obs_3, psi_obs_0, aniso_d, aniso_h, dyn_D, dyn_Q, b_odd_max, FWHM_d, sigma_d, d_fil, sigma_fil, h_cr, l_fil, r1, r2, ne_cone, ne_cone_0, h_cr_fil;


double PI, cell_x, cell_y, cell_z, D, h_B_d, h_B_h, FWHM_1, FWHM_2, alpha, beta,
    ia, iat, PA, PAt, pitch, sigma_1, sigma_2, sod, FWHM_e,
    sigma_e,
    FWHM_b,sigma_b, pi_max, E_d, E_h, E_d0, E_h0, n_CR, FWHM_h, sigma_h;


int main( void );
void setup_array ( void );
void setup_electric_field ( void );
void writeimage( char filename[], double map[207][207], double crval3, double cdelt3, double crval4 );
void writecube( char filename[], double map[207][207][207], double crval3, double cdelt3, double crval4 );
void printerror ( int status );


int main()
{

// 205
/* grid size */
    gs = 205;
    
    PI = 3.141592654;

//10 is the number for gs=200
/* cell size in arcsec */
    cell_x = 3.0;
    cell_y = 3.0;
    cell_z = 3.0;

    D = 3.94e6; //distance in pc

    h_B_d = 1600.0; //scale height disk field
    
    h_B_h = 6800.0;//scale height halo field 4.0 * 1700.0;

/* maximum of magnetic field in case of odd parity */
    
    b_odd_max = 500.0;
    
/* scaleheight of the thermal electrons */    
        
    h_e = 1400.0;// 1400.0 is from paper2, measured with halpha

    
/* electron density in cm^-3 */    

    ne_0 = 0.05;
    
    lambda_1 = 0.0637; // wavelength in meters

    lambda_2 = 0.0354; // wavelength in meters

    lambda_3 = 0.012; // wavelength in meters
        
    FWHM_1 = 6500.0; // in pc (gaussian, disk along the major axis, PaperII Fig 9 left)

    FWHM_2 = 13000.0; // in pc (gaussian, disk along the major axis)

    FWHM_h = 6500.0;// in pc (gaussian, halo along the major axis PaperII Fig 9 right)

/* FWHM of disk field in case of odd parity */    

    FWHM_d = 400; 
        
    FWHM_e = 13000.0; // electron density

/* if a vertical Gaussian magnetic field is used */     

    FWHM_b = 2000.0;
        
    sod = 5.0; // free parameter:  FWHM_1/FWHM_2 disk
        
/* disk =1: switch disk magnetic field on */

    disk = 1;
    
/* halo=1: switch halo magnetic field on */

    halo = 1;
    
/* disk_odd = 1  on */ 

    disk_odd = -1;

/* halo_odd = 1: halo field is of odd parity */

    halo_odd = 1;

/* clockwise_north = 1: if disk field is odd: disk magnetic field points in clockwise direction in the northern hemisphere */    

    clockwise_north = -1;
    

/*halo_n =1: halo mf points inward in northern hemisphere */    

    halo_n = 1;
    
/* equipartition=1: switch equipartition on */

    equipartition= -1;

/* gaussian_bfield=1: switch Gaussian vertical bfield distribution on */

    gaussian_bfield = -1;

    dynamo = -1; //?

    idyn = 1;//?
    
    dyn_D = 1.0;//?

    dyn_Q = 0.0;//?
                    
    alpha = 0.0 * PI / 180.0;// helical field

    beta = 45.0 * PI / 180.0 / 2.0;// helical field, beta_n = 45.0 * PI / 180.0;

/* pitch = pitch angle of the magnetic field spiral */

    pitch = 25.0 * PI / 180.0;

/* inclination angle; ia = 90d means edge on */    

    ia = 78.5;

    iat = ( 90.0 - ia ) * PI / 180.0;

    PA = 52.0;

    PAt = ( 90.0 - PA ) * PI / 180.0;

    E_d0 = 4.4;//4.4 is from paper2, ordered magnetic field strength disk
    
    E_h0 = 20.0;//20.0 for the filaments 4.4 is form paper2, ordered magnetic field strength halo

    aniso_d = 0.0;//0.68 is from paper2, anisotropical component
    
    aniso_h = 1.0;//0.68 is from paper2

/* conical magnetic field with filaments */

    cone = 1;

    d_fil = 40.0;//40.0

    l_fil = 300.0; //300.0

    h_cr = 800.0;

    h_cr_fil = 200.0;
    
    ne_cone_0 = 2.0;//1.0
            
    

/*************************************************************************************** Here ends the parameter input ***************/    
    
        
    setup_array ();
    setup_electric_field ();


    writeimage ( "Q0.fits", Q_0, 0.0, 0.0, 2.0 );
    writeimage ( "Q1.fits", Q_1, 4.71e9, 5.0e8, 2.0 ); // Frequenz1,Bandbreite1,?
    writeimage ( "Q2.fits", Q_2, 8.46e9, 5.0e8, 2.0 ); // Frequenz1,Bandbreite1,?
    writeimage ( "Q3.fits", Q_3, 25.0e9, 5.0e8, 2.0 );
    writeimage ( "U0.fits", U_0, 0.0, 0.0, 3.0 );
    writeimage ( "U1.fits", U_1, 4.71e9, 5.0e8, 3.0 );
    writeimage ( "U2.fits", U_2, 8.46e9, 5.0e8, 3.0 );
    writeimage ( "U3.fits", U_3, 25.0e9, 5.0e8, 3.0 );
    writeimage ( "RM.fits", RM, 4.71e9, 5.0e8, 0.0 );
    writeimage ( "RMs.fits", RM_2, 4.71e9, 5.0e8, 0.0 );


    return(0);
}

void setup_array ( void )

{

    x[0]=0;
    y[0]=0;
    z[0]=0;

    delta_x = cell_x * D * tan(PI / (3600.0 * 180.0));
    delta_y = cell_y * D * tan(PI / (3600.0 * 180.0));
    delta_z = cell_z * D * tan(PI / (3600.0 * 180.0));

    
    for (i=1; i <= gs; i++)
    {
	
	x[i] = x[i-1] + delta_x;

    }
    
    for (j=1; j <= gs; j++)
		
    {
	
	y[j]= y[j-1] + delta_y;
	
    }

    for (k=1; k <= gs; k++)
		
    {
	
	z[k]= z[k-1] + delta_z;

    }
		
    x_ref = 0.5 * x[gs];

    y_ref = 0.5 * y[gs];

    z_ref = 0.5 * z[gs];


    for (i=1; i <= gs; i++)
	
    {
	
	for (j=1; j <= gs; j++)
	
	{
	
	    for (k=1; k <= gs; k++)

	    {

	      /* x_ss[i][k] = ( x[i] - x_ref );                              A1

		y_ss[i][j][k] = cos( iat ) * ( y[j] - y_ref ) -               A2
		    sin( iat ) * ( z[k] - z_ref );

		 z_ss[i][j][k] = cos ( iat ) * ( z[k] - z_ref ) +             A3
		sin ( iat ) * ( y[j] - y_ref ); */


		x_ss[i][k] = cos ( PAt ) * ( x[i] - x_ref ) - sin( PAt ) * 
		    ( z[k] - z_ref );
		
		y_ss[i][j][k] = -sin( PAt ) * sin( iat ) * ( x[i] - x_ref ) +
		    cos( iat ) * ( y[j] - y_ref ) - sin( iat ) * 
		    cos ( PAt ) * ( z[k] - z_ref );

		z_ss[i][j][k] = sin( PAt ) * cos( iat ) * ( x[i] - x_ref ) + 
		    sin( iat ) * ( y[j] - y_ref ) +
		    cos( iat ) * cos( PAt ) * ( z[k] - z_ref );

	    }

	}

    }
    

    return;

}


void setup_electric_field ( void )

{



    if (equipartition == 1)

    {

        FWHM_1 = 2.0 * FWHM_1;
    
        FWHM_2 = 2.0 * FWHM_2;

        FWHM_b = 2.0 * FWHM_b;
    
        FWHM_h = 2.0 * FWHM_h;

        FWHM_d = 2.0 * FWHM_d;
        
    }
    
            

    sigma_1 = FWHM_1 / ( 2.0 * sqrt( 2.0 * log (2.0) ) );

    sigma_2 = FWHM_2 / ( 2.0 * sqrt( 2.0 * log (2.0) ) );

    sigma_e = FWHM_e / ( 2.0 * sqrt( 2.0 * log (2.0) ) );

    sigma_b = FWHM_b / ( 2.0 * sqrt( 2.0 * log (2.0) ) );

    sigma_h = FWHM_h / ( 2.0 * sqrt( 2.0 * log (2.0) ) );

    sigma_d = FWHM_d / ( 2.0 * sqrt( 2.0 * log (2.0) ) );

    sigma_fil = d_fil / ( 2.0 * sqrt( 2.0 * log (2.0) ) );


    for (i=1; i <= gs; i++)

    {

        for (k=1; k <= gs; k++)

        {

	    
            delta_phi_1[i][1][k] = 0.0;
            delta_phi_2[i][1][k] = 0.0;
            
//            RM[i][k] = 0.0;
            RM_2[i][k] = 0.0;
//            Q_0[i][k] = 0.0;
//            U_0[i][k] = 0.0;
            delta_phi_o1[i][k] = 0.0;


            for (j=1; j <= gs; j++)

            {

	      if (gaussian_bfield == 1) // scheibenfeldstaerke A4

                    E_d = E_d0 * exp (-pow( z_ss[i][j][k], 2.0 )  /
                                   (2.0 * pow(sigma_b, 2.0) ) ) *
                        ( exp( -(pow( x_ss[i][k], 2.0 ) +
                                 pow( y_ss[i][j][k], 2.0 ) ) / 
                             ( 2.0 * pow( sigma_1, 2.0) ) ) + sod *
                          exp( -(pow( x_ss[i][k], 2.0 ) +
                                 pow( y_ss[i][j][k], 2.0 ) ) / 
                               ( 2.0 * pow( sigma_2, 2.0) ) ) ) /
                        (1.0 + sod);
                
                else if (disk_odd == 1)

                    E_d = E_d0 * fabs( z_ss[i][j][k] ) / b_odd_max *
                        exp ( 1.0 - fabs ( z_ss[i][j][k] ) / b_odd_max );
                
                    

                    /*                   if (fabs( z_ss[i][j][k] ) <= b_odd_max )

                        E_d = E_d0 * ( pow ( sin( ( z_ss[i][j][k] /
                                                    b_odd_max ) * PI /
                                                  2.0 ), 2.0 ) );
                
                    else

                        E_d = E_d0 * exp ( -pow( fabs (z_ss[i][j][k] ) -
                                                 b_odd_max, 2.0 ) /
                                                 (2.0 * pow(sigma_b, 2.0) ) );*/
                   
                

                else
                                      
                    E_d = E_d0 * exp (-sqrt( pow( z_ss[i][j][k], 2.0 ) ) /
                                      h_B_d) * 
                        ( exp(-pow( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                          pow( y_ss[i][j][k], 2.0 ) ), 2.0 ) / 
                              ( 2.0 * pow( sigma_1, 2.0) ) ) + sod *
                          exp(-pow( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                          pow( y_ss[i][j][k], 2.0 ) ), 2.0 ) / 
                              ( 2.0 * pow( sigma_2, 2.0) ) ) ) /
                        (1.0 + sod);

               
		
                if ( halo == 1 )
                
                    E_h = E_h0 * exp (-sqrt( pow( z_ss[i][j][k], 2.0 ) ) /
                                      h_B_h) *
                        exp(-pow( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                        pow( y_ss[i][j][k], 2.0 ) ), 2.0 ) / 
                            ( 2.0 * pow( sigma_h, 2.0) ) );


                if ( dynamo == 1 )
                
                    E_h = E_h0 * exp ( - sqrt ( x_ss[i][k] * x_ss[i][k] +
                                                y_ss[i][j][k] * y_ss[i][j][k] +
                                                z_ss[i][j][k] * z_ss[i][j][k] )
                                       / h_B_h);

                if ( cone == 1 )

                {
                    
                    r1 = fabs ( z_ss[i][j][k] ) *
                        sqrt ( 1.0 / pow ( cos ( beta ), 2.0  )
                               - 1.0 )
                        + l_fil - d_fil / (2.0 * cos (beta) );

                    r2 = r1 + d_fil / cos (beta);

                    ne_cone = ne_cone_0 *
                        exp (-sqrt( pow( z_ss[i][j][k], 2.0 ) ) / h_e) *
                        exp(-pow( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                        pow( y_ss[i][j][k], 2.0 ) ), 2.0 ) / 
                            ( 2.0 * pow( sigma_e, 2.0) ) );

                    E_h = E_h * exp(-pow( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                                pow( y_ss[i][j][k], 2.0 ) )
                                          - (r1 + r2) / 2.0, 2.0 ) / 
                            ( 2.0 * pow( sigma_fil, 2.0) ) );


                    if ( fabs ( z_ss[i][j][k] ) < 1200.0 )
                        
                        alpha = 0.5 * PI / 2.0 * fabs ( z_ss[i][j][k] )
                            / 1200.0;

                    else
                        
                    alpha = 0.5 * PI / 2.0;



                    if ( sqrt( pow( x_ss[i][k], 2.0 ) + 
                               pow( y_ss[i][j][k], 2.0 ) ) < r1  )


                        ne_cone = 0.0;
                         

                    if ( sqrt( pow( x_ss[i][k], 2.0 ) + 
                               pow( y_ss[i][j][k], 2.0 ) ) > r2  )

                        ne_cone = 0.0;


                    /*                   if ( (sqrt( pow( x_ss[i][k], 2.0 ) + 
                                pow( y_ss[i][j][k], 2.0 ) ) > (r1 + r2) / 2.0 )  && fabs( z_ss[i][j][k] ) >= 400.0)

                                E_h = - E_h;   */

                
                }
                
                

                
                E_zs_h = 0.0;
                
                E_phi_h = 0.0;
                    
                E_R_h = 0.0;

                E_phi_d = 0.0;

                E_R_d = 0.0;

                E_zs_d = 0.0;


                if ( halo == 1 )

                {

                    if ( halo_odd == 1 )
                
                        E_zs_h = E_h / sqrt( 1.0 + pow( tan( alpha ), 2 ) + 
                                         pow ( tan ( beta ), 2 ) );
                    else
                
                    {
            
                        if ( ( z_ss[i][j][k] ) >=0 )
                            
                            E_zs_h = E_h / sqrt( 1.0 + pow( tan( alpha ), 2 ) + 
                                             pow ( tan ( beta ), 2 ) );
		
                        else
                    
                            E_zs_h = -E_h /
                                sqrt( 1.0 + pow( tan( alpha ), 2 ) + 
                                  pow ( tan ( beta ), 2 ) );
                        
                    }

                    if ( halo_n == 1 )

                        E_zs_h = -E_zs_h;
            
                    if ( ( z_ss[i][j][k] ) >=0 )
                        
                        E_phi_h = E_zs_h * tan( alpha );
        
                    else

                        E_phi_h = -E_zs_h * tan( alpha );
		    
                    
                    if ( ( z_ss[i][j][k] ) >=0 )

                        E_R_h = E_zs_h * tan ( beta );
                    
                    else
                
                        E_R_h = -E_zs_h * tan ( beta );

                    
                    
                }


                if ( dynamo == 1 )

                {
                    
                    theta = atan2 ( sqrt ( x_ss[i][k] * x_ss[i][k] +
                                          y_ss[i][j][k] * y_ss[i][j][k]),
                                    z_ss[i][j][k]);


                    E_R_h = 3.0 * dyn_D * sin ( theta ) * cos ( theta )
                        + 3.0 * dyn_Q * cos ( theta ) * cos ( theta ) *
                        sin ( theta ) + 3.0 * dyn_Q * sin ( theta ) *
                          ( 3.0 * cos ( theta ) * cos ( theta ) - 1.0 );

                    E_zs_h = - dyn_D * sin ( theta ) * sin ( theta ) +
                        2.0 * dyn_D * cos ( theta ) * cos ( theta )
                        - 3.0 * dyn_Q * sin ( theta ) * sin ( theta ) *
                        cos ( theta ) + 3.0 * dyn_Q * cos ( theta ) *
                        ( 3.0 * cos ( theta ) * cos ( theta ) - 1.0 );


                    E_R_h = E_h / sqrt ( E_R_h * E_R_h + E_zs_h * E_zs_h ) *
                        E_R_h;

                    E_zs_h = E_h / sqrt ( E_R_h * E_R_h + E_zs_h * E_zs_h ) *
                    E_zs_h;

                    if ( idyn == 1 )

                    {

                        E_R_h = - E_R_h;

                        E_zs_h = - E_zs_h;
                        
                    }
                    


                }
                


                if ( disk == 1) 

                {

		    
                    if ( disk_odd == 1 )

                    {
			
                        if ( ( z_ss[i][j][k] ) >=0 ) 

                        {
			
                            E_phi_d = E_d / sqrt( 1.0 + pow( tan( pitch ), 2 ) );
                
                        }
			
                        else
                
                        {

                            E_phi_d = -E_d / sqrt( 1.0 + pow( tan( pitch ), 2 ) );

                        }

                    }

                    else
                
                    {

                        E_phi_d = E_d / sqrt( 1.0 + pow( tan( pitch ), 2 ) );

                   }

                   
                    E_R_d = - tan( pitch ) * E_phi_d;
                    
                }

                
                if ( clockwise_north == 1 )

                {
                    
                    E_phi_d = -E_phi_d;

                    E_R_d = -E_R_d;
                        
                }
 

                
                phi = atan2( y_ss[i][j][k], x_ss[i][k] );



                E_xs_d = E_R_d * cos( phi ) - E_phi_d * sin( phi );
		
                E_ys_d = E_R_d * sin( phi ) + E_phi_d * cos( phi );

                E_x_d[i][j][k] = cos( PAt ) * E_xs_d - sin( PAt ) *
                    sin( iat ) * E_ys_d + sin( PAt ) * cos ( iat ) * E_zs_d;

                E_y_d[i][j][k] = cos( iat ) * E_ys_d + sin( iat ) * E_zs_d;
                
                E_z_d[i][j][k] = -sin( PAt ) * E_xs_d - sin( iat ) *
                    cos( PAt ) * E_ys_d + cos( iat ) * cos( PAt ) * E_zs_d;


                
                E_xs_h = E_R_h * cos( phi ) - E_phi_h * sin( phi );
		
                E_ys_h = E_R_h * sin( phi ) + E_phi_h * cos( phi );

                E_x_h[i][j][k] = cos( PAt ) * E_xs_h - sin( PAt ) *
                    sin( iat ) * E_ys_h + sin( PAt ) * cos ( iat ) * E_zs_h;

                E_y_h[i][j][k] = cos( iat ) * E_ys_h + sin( iat ) * E_zs_h;
                
                E_z_h[i][j][k] = -sin( PAt ) * E_xs_h - sin( iat ) *
                    cos( PAt ) * E_ys_h + cos( iat ) * cos( PAt ) * E_zs_h;

               
                psi_d[i][j][k] = atan2(-E_x_d[i][j][k], E_z_d[i][j][k]) -
                    PI / 2.0;

                psi_h[i][j][k] = atan2(-E_x_h[i][j][k], E_z_h[i][j][k]) -
                    PI / 2.0;

                ne = ne_0 * exp (-sqrt( pow( z_ss[i][j][k], 2.0 ) ) / h_e) *
                        exp(-pow( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                        pow( y_ss[i][j][k], 2.0 ) ), 2.0 ) / 
                            ( 2.0 * pow( sigma_e, 2.0) ) );
                
                delta_phi_1[i][j+1][k] = delta_phi_1[i][j][k] -
                        0.81 * pow( lambda_1, 2.0) *
                    ( ne + ne_cone ) * (aniso_d * E_y_d[i][j][k] + aniso_h *
                    E_y_h[i][j][k]) * delta_y;


 
                delta_phi_2[i][j+1][k] = delta_phi_2[i][j][k] - 
                    0.81 * pow( lambda_2, 2.0) *
                    ( ne + ne_cone ) * (aniso_d * E_y_d[i][j][k] + aniso_h *
                     E_y_h[i][j][k]) * delta_y;

                delta_phi_3[i][j+1][k] = delta_phi_3[i][j][k] - 
                    0.81 * pow( lambda_3, 2.0) *
                    ( ne + ne_cone ) * (aniso_d * E_y_d[i][j][k] + aniso_h *
                     E_y_h[i][j][k]) * delta_y;


                delta_phi_o1[i][k] = delta_phi_1[i][10][k] -
                    delta_phi_2[i][10][k];


/* RM for a Faraday Screen. The "-" corrects that the observer is at the front side of the cube. */

                RM_2[i][k] = RM_2[i][k] - 0.81 *
                    ( ne + ne_cone ) * (aniso_d * E_y_d[i][j][k] + aniso_h *
                     E_y_h[i][j][k])* delta_y;

  

                                
            }

        }
    
    }



    for (i=1; i <= gs; i++)

    {

        for (k=1; k <= gs; k++)
	    
        {

            Q_0[i][k] = 0.0;
	    
            U_0[i][k] = 0.0;
	    
            Q_1[i][k] = 0.0;

            U_1[i][k] = 0.0;
	    
            Q_2[i][k] = 0.0;

            U_2[i][k] = 0.0;

            Q_3[i][k] = 0.0;

            U_3[i][k] = 0.0;

            Et[i][k] = 0.0;

            RM[i][k] = 0.0;

            
            for (j=1; j <= gs; j++)
                
            {

                
                psi_1_d = psi_d[i][j][k] + delta_phi_1[i][j][k];
		
                psi_2_d = psi_d[i][j][k] + delta_phi_2[i][j][k];

                psi_3_d = psi_d[i][j][k] + delta_phi_3[i][j][k];

                psi_1_h = psi_h[i][j][k] + delta_phi_1[i][j][k];
		
                psi_2_h = psi_h[i][j][k] + delta_phi_2[i][j][k];

                psi_3_h = psi_h[i][j][k] + delta_phi_3[i][j][k];

                 if (equipartition == 1)

                     n_CR = (pow( E_x_d[i][j][k], 2.0 ) +
                             pow( E_y_d[i][j][k], 2.0 ) +
                             pow( E_z_d[i][j][k], 2.0 ) +
                             pow( E_x_h[i][j][k], 2.0 ) +
                             pow( E_y_h[i][j][k], 2.0 ) +
                             pow( E_z_h[i][j][k], 2.0 ) ) /
                         pow(E_d0 + E_h0, 2.0);

                 else if ( cone == 1 )

                {
                    
                    r1 = fabs ( z_ss[i][j][k] ) *
                        sqrt ( 1.0 / pow ( cos ( beta ), 2.0  )
                               - 1.0 )
                        + l_fil - d_fil / (2.0 * cos (beta) );

                    r2 = r1 + d_fil / cos (beta);

                    if ( ( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                 pow( y_ss[i][j][k], 2.0 ) ) >= r1 )
                        && ( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                   pow( y_ss[i][j][k], 2.0 ) ) <= r2 ) )

                        n_CR =  exp (-sqrt( pow( z_ss[i][j][k], 2.0 ) ) /
                                     h_cr_fil) *
                            ( exp(-pow( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                              pow( y_ss[i][j][k], 2.0 ) ),
                                        2.0 ) / 
                                  ( 2.0 * pow( sqrt ( 2.0 ) *
                                               sigma_1, 2.0) ) ) + sod *
                              exp(-pow( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                              pow( y_ss[i][j][k], 2.0 ) ),
                                        2.0 ) / 
                                  ( 2.0 * pow( sqrt ( 2.0 ) * sigma_2,
                                               2.0) ) ) ) / (1.0 + sod);
                    else

                        n_CR =  exp (-sqrt( pow( z_ss[i][j][k], 2.0 ) ) /
                                     h_cr) * 
                            ( exp(-pow( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                              pow( y_ss[i][j][k], 2.0 ) ),
                                        2.0 ) / ( 2.0 * pow( sqrt ( 2.0 ) *
                                               sigma_1, 2.0) ) ) + sod *
                              exp(-pow( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                              pow( y_ss[i][j][k], 2.0 ) ),
                                        2.0 ) / 
                                  ( 2.0 * pow( sqrt ( 2.0 ) * sigma_2,
                                               2.0) ) ) ) / (1.0 + sod);
                         
                }
                 
                 else

                     n_CR =  exp (-sqrt( pow( z_ss[i][j][k], 2.0 ) ) /
                                  h_cr) * 
                         ( exp(-pow( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                           pow( y_ss[i][j][k], 2.0 ) ), 2.0 ) / 
                               ( 2.0 * pow( sqrt ( 2.0 ) *
                                            sigma_1, 2.0) ) ) + sod *
                           exp(-pow( sqrt( pow( x_ss[i][k], 2.0 ) + 
                                           pow( y_ss[i][j][k], 2.0 ) ), 2.0 ) / 
                               ( 2.0 * pow( sqrt ( 2.0 ) * sigma_2, 2.0) ) ) ) /
                         (1.0 + sod);
                 

                Q_0[i][k]= Q_0[i][k] + n_CR *
                    ( ( pow( E_x_d[i][j][k], 2.0) +
                        pow( E_z_d[i][j][k], 2.0) ) *
                      cos(2.0 * psi_d[i][j][k]) +
                      ( pow( E_x_h[i][j][k], 2.0) +
                        pow( E_z_h[i][j][k], 2.0) ) *
                      cos(2.0 * psi_h[i][j][k]) );

                U_0[i][k]= U_0[i][k] + n_CR *
                    ( ( pow( E_x_d[i][j][k], 2.0) +
                        pow( E_z_d[i][j][k], 2.0) ) *
                      sin(2.0 * psi_d[i][j][k]) +
                      ( pow( E_x_h[i][j][k], 2.0) +
                        pow( E_z_h[i][j][k], 2.0) ) *
                      sin(2.0 * psi_h[i][j][k]) );

 
                Q_1[i][k]= Q_1[i][k] + n_CR *
                    ( ( pow( E_x_d[i][j][k], 2.0) +
                        pow( E_z_d[i][j][k], 2.0) ) *
                      cos(2.0 * psi_1_d) +
                      ( pow( E_x_h[i][j][k], 2.0) +
                        pow( E_z_h[i][j][k], 2.0) ) *
                        cos(2.0 * psi_1_h) );

                U_1[i][k]= U_1[i][k] + n_CR *
                    ( ( pow( E_x_d[i][j][k], 2.0) +
                        pow( E_z_d[i][j][k], 2.0) ) *
                      sin(2.0 * psi_1_d) +
                      ( pow( E_x_h[i][j][k], 2.0) +
                        pow( E_z_h[i][j][k], 2.0) ) *
                        sin(2.0 * psi_1_h) );


                Q_2[i][k]= Q_2[i][k] + n_CR *
                    ( ( pow( E_x_d[i][j][k], 2.0) +
                        pow( E_z_d[i][j][k], 2.0) ) *
                      cos(2.0 * psi_2_d) +
                      ( pow( E_x_h[i][j][k], 2.0) +
                        pow( E_z_h[i][j][k], 2.0) ) *
                        cos(2.0 * psi_2_h) );


                U_2[i][k]= U_2[i][k] + n_CR *
                    ( ( pow( E_x_d[i][j][k], 2.0) +
                        pow( E_z_d[i][j][k], 2.0) ) *
                      sin(2.0 * psi_2_d) +
                      ( pow( E_x_h[i][j][k], 2.0) +
                        pow( E_z_h[i][j][k], 2.0) ) *
                      sin(2.0 * psi_2_h) );

                Q_3[i][k]= Q_3[i][k] + n_CR *
                    ( ( pow( E_x_d[i][j][k], 2.0) +
                        pow( E_z_d[i][j][k], 2.0) ) *
                      cos(2.0 * psi_3_d) +
                      ( pow( E_x_h[i][j][k], 2.0) +
                        pow( E_z_h[i][j][k], 2.0) ) *
                      cos(2.0 * psi_3_h) );

                U_3[i][k]= U_3[i][k] + n_CR *
                    ( ( pow( E_x_d[i][j][k], 2.0) +
                        pow( E_z_d[i][j][k], 2.0) ) *
                      sin(2.0 * psi_3_d) +
                      ( pow( E_x_h[i][j][k], 2.0) +
                        pow( E_z_h[i][j][k], 2.0) ) *
                      sin(2.0 * psi_3_h) );

                
                
                
                Et[i][k] = Et[i][k] + n_CR *
                    (pow( E_x_d[i][j][k], 2.0 ) + pow( E_y_d[i][j][k], 2.0 ) +
                     pow( E_z_d[i][j][k], 2.0 ) + pow( E_x_h[i][j][k], 2.0 ) +
                     pow( E_y_h[i][j][k], 2.0 ) + pow( E_z_h[i][j][k], 2.0 ) );

                
                RM[i][k] = 2.0 * psi_1_h;
                
                
                    
            }
	
        }
        
    }

    /* Search for maximum of polarrized intensity */

    pi_max = 0.0;
    
    for (i=1; i <= gs; i++)
        
    {

        for (k=1; k <= gs; k++)
            
        {

            PI_0[i][k] = sqrt( pow( Q_0[i][k], 2.0) + pow( U_0[i][k], 2.0) );

            if (pi_max < PI_0[i][k])

                pi_max = PI_0[i][k];
            
        
        }


    }
   


    for (i=1; i <= gs; i++)

    {

        for (k=1; k <= gs; k++)
            
        {

            psi_obs_0 = 0.5 * atan2( U_0[i][k], Q_0[i][k] );

/* use this when you want to use the line-of-sight integrated RM */

/*            psi_obs_1 = psi_obs_0 + RM[i][k] * pow(lambda_1, 2.0);

            psi_obs_2 = psi_obs_0 + RM[i][k] * pow(lambda_2, 2.0);
            
            Q_1[i][k] = PI_0[i][k] * cos(2.0 * psi_obs_1);

            U_1[i][k] = PI_0[i][k] * sin(2.0 * psi_obs_1);

            Q_2[i][k] = PI_0[i][k] * cos(2.0 * psi_obs_2);

            U_2[i][k] = PI_0[i][k] * sin(2.0 * psi_obs_2);*/

/* compute RM internal in the code instead of using MIRIAD/AIPS */

            psi_obs_1 = 0.5 * atan2( U_1[i][k], Q_1[i][k] );

            psi_obs_2 = 0.5 * atan2( U_2[i][k], Q_2[i][k] );

            
            PI_1[i][k] = sqrt( pow( Q_1[i][k], 2.0) + pow( U_1[i][k], 2.0) );
            
            PI_2[i][k] = sqrt( pow( Q_2[i][k], 2.0) + pow( U_2[i][k], 2.0) );
            



        }

    }


    for (i=1; i <= gs; i++)

    {

        for (k=1; k <= gs; k++)
	    
        {

            Q_0[i][k] = 1.0e-6 * Q_0[i][k];
	    
            U_0[i][k] = 1.0e-6 * U_0[i][k];
	    
            Q_1[i][k] = 1.0e-6 * Q_1[i][k];

            U_1[i][k] = 1.0e-6 * U_1[i][k];
	    
            Q_2[i][k] = 1.0e-6 * Q_2[i][k];

            U_2[i][k] = 1.0e-6 * U_2[i][k];

            Q_3[i][k] = 1.0e-6 * Q_3[i][k];

            U_3[i][k] = 1.0e-6 * U_3[i][k];


        }


    }
    


  

    return;

}

void writeimage(char filename[], double map[207][207], double crval3, double cdelt3, double crval4 )

    /******************************************************/
    /* Create a FITS primary array containing a 2-D image */
    /******************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii, jj;
    long  nelements;
    long fpixel;
    double bscale, bzero, epoch, obsra, obsdec, crval1, cdelt1, crpix1, crota1,
        crval2, cdelt2, crpix2, crota2, crpix3, crota3, cdelt4, crpix4;
    int blank;
    double *array[206];

    /* initialize FITS image parameters */
/*    char filename[] = "Q.fits"; name for new FITS file */
    int bitpix   =  FLOAT_IMG; /* 16-bit unsigned short pixel values       */
    long naxis    =   4;  /* 4-dimensional image                           */
    long naxes[4] = { 205, 205, 1, 1 };   /* image is 300 pixels wide by 200 rows */

    /* allocate memory for the whole image */
    array[0] = (double *)malloc( naxes[0] * naxes[1]
                                        * sizeof( double ) );

//    /* initialize pointers to the start of each row of the image */
    for( ii=1; ii<naxes[1]; ii++ )
      array[ii] = array[ii-1] + naxes[0];

    remove(filename);               /* Delete old file if it already exists */

    status = 0;         /* initialize status before calling fitsio routines */

    if (fits_create_file(&fptr, filename, &status)) /* create new FITS file */
         printerror( status );           /* call printerror if error occurs */

    /* write the required keywords for the primary array image.     */
    /* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = 16 (signed short integers) with   */
    /* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
    /* FITS uses to store unsigned integers.  Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */

    if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
         printerror( status );

    /* initialize the values in the image with a linear ramp function */

    for (jj = 0; jj < naxes[1]; jj++)
    {   for (ii = 0; ii < naxes[0]; ii++)
        {
            array[jj][ii] = (ii+jj) - (ii+jj) + map[ii][jj];
        }
    }

    fpixel = 1;                               /* first pixel to write      */
    nelements = naxes[0] * naxes[1];          /* number of pixels to write */

    /* write the array of unsigned integers to the FITS file */
    if ( fits_write_img(fptr, TDOUBLE, fpixel, nelements, array[0], &status) )
        printerror( status );
      
    free( array[0] );  /* free previously allocated memory */

    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
/*    naxis = 3;
    if ( fits_update_key(fptr, TINT, "NAXIS", &naxis,
			 "number of data axis", &status) )
			 printerror( status ); */

    if ( fits_update_key(fptr, TSTRING, "OBJECT", "NGC 253", "Source name",
			 &status) )
         printerror( status );

    if ( fits_update_key(fptr, TSTRING, "TELESCOPE", "Model", NULL,
			 &status) )
         printerror( status );

    bscale = 1.0;
    if ( fits_update_key(fptr, TDOUBLE, "BSCALE", &bscale,
			 "REAL = TAPE * BSCALE + BZERO", &status) )
         printerror( status );

    bzero = 0.0;
    if ( fits_update_key(fptr, TDOUBLE, "BZERO", &bzero, NULL, &status) )
         printerror( status );

    if ( fits_update_key(fptr, TSTRING, "BUNIT", "tbd", "Units of flux",
			 &status) )
         printerror( status );

    epoch = 2.0e3;
    if ( fits_update_key(fptr, TDOUBLE, "EPOCH", &epoch, NULL,
			 &status) )
         printerror( status );

    blank = -1;
    if ( fits_update_key(fptr, TINT, "BLANK", &blank,
			 "IEEE not-a-number for blanked pixels", &status) )
         printerror( status );

    obsra = 1.18870383535e1;
    if ( fits_update_key(fptr, TDOUBLE, "OBSRA", &obsra, "Antenna pointing RA",
			 &status) )
         printerror( status );

    obsdec = -2.52886438077e1;
    if ( fits_update_key(fptr, TDOUBLE, "OBSDEC", &obsdec,
			 "Antenna pointing DEC", &status) )
	printerror( status );


    if ( fits_update_key(fptr, TSTRING, "CTYPE1", "RA---SIN", NULL,
			 &status) )
         printerror( status );

    crval1 =  1.18880416667e1;
    if ( fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1, NULL,
			 &status) )
         printerror( status );

    cdelt1 =  -2.777777845e-3 * cell_x / 10.0;
    if ( fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1, NULL,
			 &status) )
         printerror( status );

    crpix1 =  1.031999969e2;
    if ( fits_update_key(fptr, TDOUBLE, "CRPIX1", &crpix1, NULL,
			 &status) )
         printerror( status );

    crota1 =  0.0;
    if ( fits_update_key(fptr, TDOUBLE, "CROTA1", &crota1, NULL,
			 &status) )
         printerror( status );

    if ( fits_update_key(fptr, TSTRING, "CTYPE2", "DEC--SIN", NULL,
			 &status) )
         printerror( status );

    crval2 =  -2.52882777778e1;
    if ( fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2, NULL,
			 &status) )
         printerror( status );

    cdelt2 =  2.777777845e-3 * cell_y / 10.0;
    if ( fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2, NULL,
			 &status) )
         printerror( status );

    crpix2 =  1.031999969e2;
    if ( fits_update_key(fptr, TDOUBLE, "CRPIX2", &crpix2, NULL,
			 &status) )
         printerror( status );

    crota2 =  0.0;
    if ( fits_update_key(fptr, TDOUBLE, "CROTA2", &crota2, NULL,
			 &status) )
         printerror( status );

    if ( fits_update_key(fptr, TSTRING, "CTYPE3", "FREQ", NULL,
			 &status) )
         printerror( status );

/*    crval3 =  4.85e9; */
    if ( fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3, NULL,
			 &status) )
         printerror( status );

/*    cdelt3 =  5.0e8; */
    if ( fits_update_key(fptr, TDOUBLE, "CRDELT3", &cdelt3, NULL,
			 &status) )
         printerror( status );

    crpix3 =  1.0;
    if ( fits_update_key(fptr, TDOUBLE, "CRPIX3", &crpix3, NULL,
			 &status) )
         printerror( status );

    crota3 =  0.0;
    if ( fits_update_key(fptr, TDOUBLE, "CROTA3", &crota3, NULL,
			 &status) )
         printerror( status );

    if ( fits_update_key(fptr, TSTRING, "CTYPE4", "STOKES", NULL,
                         &status) )
        printerror( status );

    /* crval4 = 2.0: Stokes Q, crval4 = 3.0: Stokes U */
    if ( fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4, NULL,
			 &status) )
         printerror( status );

    cdelt4 =  1.0;
    if ( fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4, NULL,
			 &status) )
         printerror( status );

    crpix4 =  1.0;
    if ( fits_update_key(fptr, TDOUBLE, "CRPIX4", &crpix4, NULL,
			 &status) )
         printerror( status );


    if ( fits_close_file(fptr, &status) )                /* close the file */
         printerror( status );

    return;
}


void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/


    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}
