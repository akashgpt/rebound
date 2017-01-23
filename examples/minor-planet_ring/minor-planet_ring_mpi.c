/**
 * Minor-planet ring
 *
 * Simulatin a ring/disk of particles surroinding a triaxial ellipsoidal
 * body (irregular later) in a reference frame rotating with the central
 * body (uniform - straight forward to extend to arbitrary rotation of 
 * central body).
 * 
 */
 
 
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include <gsl/gsl_integration.h>
#include <time.h>


const double G = 6.67384E-11;//Gravitational constant

//Force due to a triaxial ellipsoidal at an external point (x,y,z)
//: Fx=beta*Bx*x; Fy=beta*By*y; Fz=beta*Bz*z. 
double Bx;
double By;
double Bz;
double V_ell;//Potential of triaxial ellipsoid


void cb_force(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);

//function call to calculate Bx, By, Bz
int int_ellipsoid (double x, double y, double z);

double coefficient_of_restitution(const struct reb_simulation*r, double v){
     return 0.5;
}

double tmax = 1E4;


//central-body (cb) parameters (>> if want to remove, add as reqd input 
//argument for int_ellipsoid, except omega)
const double sma_ax=124.E3;//semi-major axis in m (ax>=ay>=az)
const double sma_ay=124.E3;
const double sma_az=124.E3;
const double rho_cb=1.840E3;//cb densiy in kg/m3
const double rot_period = 25200.;//cd rotation period in secs (e.g.: ~7hrs for Chariklo)

const double rho_ptcl = 900.;//particle density > assuming single





// MAIN $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
int main(int argc, char* argv[]){

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//Files written during program

	FILE *fpos = fopen("pos.txt","w");
	FILE *fvel = fopen("vel.txt","w");
	FILE *fprop = fopen("prop.txt","w");
	FILE *facc = fopen("acc.txt","w");
	
	fprintf(fpos, "Ptcl-ID\t Time\t r\t rx\t ry\t rz\n");
	fprintf(fpos, "\n");
	fprintf(fvel, "Ptcl-ID\t Time\t v\t vx\t vy\t vz\n");
	fprintf(fvel, "\n");
	fprintf(facc, "Ptcl-ID\t Time\t a\t ax\t ay\t az\n");
	fprintf(facc, "\n");
	fprintf(fprop, "Time\t L\t PE\t KE\t C\n");
	fprintf(fprop, "\n");

 	fclose(fpos);
    fclose(fvel);
    fclose(facc);
    fclose(fprop);


    
    // Calculate the time taken ...
    clock_t t;
    t = clock();
    printf("Simulation started ...\n");
    

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//Simulation parameters and initialization
	struct reb_simulation* r = reb_create_simulation();
	// setup constants
	//~ tmax							=	1.*rot_period;	
	r->G 							= 6.67384E-11;
	r->dt 							= 1.E-10; // initial timestep
	r->integrator					= REB_INTEGRATOR_IAS15; //>>HERMES?
	//~ r->integrator					= REB_INTEGRATOR_LEAPFROG; //>>HERMES?
	r->ri_ias15.epsilon 			= 1.E-4; // accuracy parameter
	r->gravity						= REB_GRAVITY_TREE;
	r->opening_angle2				= 1.5; // This constant determines the accuracy of the tree code gravity estimate.
	//~ r->boundary						= REB_BOUNDARY_OPEN;
	const double boxsize 			= 1000.E3;
	reb_configure_box(r,boxsize,2,2,1);
	//~ r->softening 	= 0.02;		// Gravitational softening length
	r->collision 					= REB_COLLISION_TREE;
	r->coefficient_of_restitution   = coefficient_of_restitution;
	r->force_is_velocity_dependent	= 1;
	r->additional_forces			= cb_force;	// setup callback function for velocity dependent forces
	r->heartbeat					= heartbeat;
	
	
    // Initialize MPI
    // This can only be done after reb_configure_box.
    reb_mpi_init(r);

	
	//If want to run for a single particle, comment the particle 
	//initialization loops below and uncomment the immediate statements
	//>>>
	//~ struct reb_particle pt = {0};
	//~ pt.x = 10000000;
	//~ pt.x = 4.224162093953768e7;
	//~ pt.vy = 5.586378386222478e3;
	//~ pt.vy = -4.547473508864641e-13;//5.586378386222478e3;
	//~ pt.m = 1.;
	//~ reb_add(r, pt);
	//<<<


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	// Particle disk parameters
	
	//total number of particles in ring/disk: MULTIPLE OF 8 required 
	//under the current scheme, will put checks later to satisfy this 
	//irrespctive of the input
    int N_ring_temp = 800;
    if(N_ring_temp%8 != 0){
		N_ring_temp = -N_ring_temp%8 + N_ring_temp + 8; 
	} 
	while(N_ring_temp%r->mpi_num != 0){
		N_ring_temp = N_ring_temp + 8;
	}
	printf("Total ring particles= %d\n",N_ring_temp);
	int	N_ring = N_ring_temp/r->mpi_num;
	//~ printf("N_ring per core= %d\n",N_ring);
    r->N_active = N_ring;
      
    double a_min = 200.E3, a_max = 400.E3; //radial extent of ring
    double theta_min = 0., theta_max = M_PI/2.; //for uniform initialization of particles
    double disk_thickness = 10.E3;
    double ptcl_rad = 1.e3; //particle radius - assuming const. for now
    double ptcl_mass = 0.*((4./3.)*M_PI*sma_ax*sma_ay*sma_az*rho_ptcl);//particle mass
	
	//initialized velocity is circular, Keplerian assuming a spherical 
	//body at (0,0,0) with the mass of the actual cb; velfld_param is 
	//there to tweak and experiment with this velocity field
	//SUGGESTIONS ? 
    double velfld_param = 1.0;
    double eff_cb_mass = velfld_param*((4./3.)*M_PI*sma_ax*sma_ay*sma_az*rho_cb);

	double *x_ptr, *y_ptr, *z_ptr;
    x_ptr = (double*) calloc(N_ring, sizeof(double));
    y_ptr = (double*) calloc(N_ring, sizeof(double));
    z_ptr = (double*) calloc(N_ring, sizeof(double));

    double powerlaw = 1;
    int i_ptcl = 0;
    int j_ptcl = 0;
    int inserted = 0;
    FILE *fptcl_init = fopen("ptcl_initialization.txt","w+");

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	//particle initialization
    while(r->N < N_ring){
		i_ptcl=i_ptcl+1;
		inserted = 0;
		double x_ptcl, y_ptcl, z_ptcl, a_ptcl, theta_ptcl;
		
        while(inserted == 0){
	        a_ptcl = reb_random_powerlaw(a_min, a_max, powerlaw);
	        theta_ptcl = reb_random_uniform(theta_min, theta_max);
	        x_ptcl = a_ptcl*cos(theta_ptcl);
	        y_ptcl = a_ptcl*sin(theta_ptcl);
	        z_ptcl = reb_random_uniform(0., disk_thickness);
	        *(x_ptr + i_ptcl-1) = x_ptcl;
	        *(y_ptr + i_ptcl-1) = y_ptcl;
	        *(z_ptr + i_ptcl-1) = z_ptcl;
	        fprintf(fptcl_init,"i_ptcl-%d -> a: %f; theta: %f; x: %f;  y: %f;  z: %f.\n",i_ptcl, a_ptcl/1000., theta_ptcl*(180./M_PI), x_ptcl/1000., y_ptcl/1000., z_ptcl/1000.);
	        
	        inserted =1;
	        
	        //Checking for overlapp with its image
			if(z_ptcl < ptcl_rad || y_ptcl < ptcl_rad || x_ptcl < ptcl_rad){
					inserted =0;
			}

	        //Checking for overlapp with other particles: 
	        //DOES your initialization scheme ensure that this is taken care of?	        
			if(i_ptcl>1 && inserted==1){
				j_ptcl=1;
				while(inserted ==1 && j_ptcl<i_ptcl){
					double dz = fabs(z_ptcl - *(z_ptr + j_ptcl-1));
					if(dz < 2.*ptcl_rad){
						double dy = fabs(y_ptcl - *(y_ptr + j_ptcl-1));
						if(dy < 2.*ptcl_rad){
							double dx = fabs(x_ptcl - *(x_ptr + j_ptcl-1));
							if(dx < 2.*ptcl_rad){
								if(sqrt(dz*dz + dy*dy + dx*dx) < 2.*ptcl_rad){
									inserted =0;
								}
							}
						}
					}
					j_ptcl=j_ptcl+1;		
				}
			}
			
			
		}
    
		double vel = sqrt(eff_cb_mass*G)/(x_ptcl*x_ptcl + y_ptcl*y_ptcl + z_ptcl*z_ptcl);
        double vx_ptcl = -vel*sin(theta_ptcl); double vy_ptcl = vel*cos(theta_ptcl); double vz_ptcl=0.;
  
        //EIGHT copies of each randomized position owing to symmetry of triaxial bodies
        //particle 1 {+x +y +z}
        struct reb_particle pt1 = {0};
        pt1.x = x_ptcl; pt1.y = y_ptcl; pt1.z = z_ptcl; 
        pt1.vx = vx_ptcl; pt1.vy = vy_ptcl; pt1.vz = vz_ptcl; 
        pt1.r = ptcl_rad;
        pt1.m  = ptcl_mass;
		reb_add(r, pt1);
        //particle 2 {+x +y -z}
        struct reb_particle pt2 = {0};
        pt2.x = x_ptcl; pt2.y = y_ptcl; pt2.z = -z_ptcl; 
        pt2.vx = vx_ptcl; pt2.vy = vy_ptcl; pt2.vz = -vz_ptcl; 
        pt2.r = ptcl_rad;
        pt2.m  = ptcl_mass;
		reb_add(r, pt2);
        //particle 3 {-x +y +z}
        struct reb_particle pt3 = {0};
        pt3.x = -x_ptcl; pt3.y = y_ptcl; pt3.z = z_ptcl; 
        pt3.vx = vx_ptcl; pt3.vy = -vy_ptcl; pt3.vz = vz_ptcl; 
        pt3.r = ptcl_rad;
        pt3.m  = ptcl_mass;
		reb_add(r, pt3);
        //particle 4 {-x +y -z}
        struct reb_particle pt4 = {0};
        pt4.x = -x_ptcl; pt4.y = y_ptcl; pt4.z = -z_ptcl; 
        pt4.vx = vx_ptcl; pt4.vy = -vy_ptcl; pt4.vz = -vz_ptcl; 
        pt4.r = ptcl_rad;
        pt4.m  = ptcl_mass;
		reb_add(r, pt4);
        //particle 5 {+x -y +z}
        struct reb_particle pt5 = {0};
        pt5.x = x_ptcl; pt5.y = -y_ptcl; pt5.z = z_ptcl; 
        pt5.vx = -vx_ptcl; pt5.vy = vy_ptcl; pt5.vz = vz_ptcl; 
        pt5.r = ptcl_rad;
        pt5.m  = ptcl_mass;
		reb_add(r, pt5);
        //particle 6 {+x -y -z}
        struct reb_particle pt6 = {0};
        pt6.x = x_ptcl; pt6.y = -y_ptcl; pt6.z = -z_ptcl; 
        pt6.vx = -vx_ptcl; pt6.vy = vy_ptcl; pt6.vz = -vz_ptcl; 
        pt6.r = ptcl_rad;
        pt6.m  = ptcl_mass;
		reb_add(r, pt6);
        //particle 7 {-x -y +z}
        struct reb_particle pt7 = {0};
        pt7.x = -x_ptcl; pt7.y = -y_ptcl; pt7.z = z_ptcl; 
        pt7.vx = -vx_ptcl; pt7.vy = -vy_ptcl; pt7.vz = vz_ptcl; 
        pt7.r = ptcl_rad;
        pt7.m  = ptcl_mass;
		reb_add(r, pt7);
        //particle 8 {-x -y -z}
        struct reb_particle pt8 = {0};
        pt8.x = -x_ptcl; pt8.y = -y_ptcl; pt8.z = -z_ptcl; 
        pt8.vx = -vx_ptcl; pt8.vy = -vy_ptcl; pt8.vz = -vz_ptcl; 
        pt8.r = ptcl_rad;
        pt8.m  = ptcl_mass;
		reb_add(r, pt8);
	
	}
	
	fclose(fptcl_init);
	free(x_ptr); free(y_ptr); free(z_ptr);

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



#ifdef OPENGL
    // Hack to artificially increase particle array.
    // This cannot be done once OpenGL is activated. 
    r->allocatedN *=8;
    r->particles = realloc(r->particles,sizeof(struct reb_particle)*r->allocatedN);
#endif // OPENGL
    //~ printf("N_active= %d\n",r->N_active);
	// Start the integration
    reb_integrate(r, tmax);

    // Cleanup
    reb_mpi_finalize(r);
	reb_free_simulation(r); 

	//~ system("rm -v radius.txt");	// remove previous output
	
	t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Simulation finished. Time taken to execute: %f s\n", time_taken);
	
}
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$












//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//Central body force function
void cb_force(struct reb_simulation* r){
	
	double omega_cb=(2*M_PI/rot_period);
	double beta_cb = (2*M_PI*rho_cb*G);
	struct reb_particle* particles = r->particles;
	const int N = r->N;
	//~ const struct reb_particle star = particles[0];				// cache
#pragma omp parallel for
	for (int i=0;i<N;i++){
		int_ellipsoid(particles[i].x, particles[i].y, particles[i].z);

		double fx_cb = -beta_cb*particles[i].x*Bx;
		double fy_cb = -beta_cb*particles[i].y*By;
		double fz_cb = -beta_cb*particles[i].z*Bz;
		
		//Equations of motion
		particles[i].ax = particles[i].ax + fx_cb + (omega_cb*particles[i].vy) + (omega_cb*omega_cb*particles[i].x);
		particles[i].ay = particles[i].ay + fy_cb - (omega_cb*particles[i].vx) + (omega_cb*omega_cb*particles[i].y);
		particles[i].az = particles[i].az + fz_cb;
			
		//~ fprintf(facc, "QQQ %d\t%.15e\t%.15f\t%.15f\t%.15f\n",i+1,r->t, particles[i].ax,particles[i].ay, particles[i].az);  
		//~ fprintf(facc, "qqq %d\t%.15e\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n",i+1,r->t, fx_cb, fy_cb, fz_cb, (omega_cb*omega_cb*particles[i].x),(omega_cb*omega_cb*particles[i].y) );    
  
		//~ fclose(facc);

	}
}


// Better would be to have the below part in another file//

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//Functional form of Bx, By, Bz (ref: S. Charndrashekar's book, Ellipsoidal figures of equilibrium)
struct my_f_params { double AX0; double AY0; double AZ0; double X_pow; double Y_pow; double Z_pow;};

double my_f (double u, void * p) {
	
	struct my_f_params * params = (struct my_f_params *)p;
	double AX0 = (params->AX0);
	double AY0 = (params->AY0);
	double AZ0 = (params->AZ0);
	double X_pow = (params->X_pow);
	double Y_pow = (params->Y_pow);
	double Z_pow = (params->Z_pow);

	return ( AX0*AY0*AZ0/( pow((AX0*AX0 + u), (0.5+X_pow))*pow((AY0*AY0 + u), (0.5+Y_pow))*pow((AZ0*AZ0 + u), (0.5+Z_pow)) ) );

}

//Relative error and absolute
const double epsabs0 = 1.E-10;
const double epsrel0 = 1.E-10;
	
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//function call to calculate Bx, By, Bz: utilizes GSL's functions
int int_ellipsoid (double x, double y, double z) {
	
	//non-dimesionalized parameters: to avoid numerical issues with operations on large numbers
	double x0, y0, z0;
	double ax0, ay0, az0;
	
	double a_nd;

	double p, q, r;  //cubic eqn's coeff.: solving for limit
	double	sol_cubic[3];
	double	lambda=0.;
	double result[3];
	double error[3];
	double pow_xyz[3] = {0.,0.,0.};
	size_t limit0=1000000;
	int i;
	
	a_nd=sma_ax;
	
	x0=x/a_nd;
	y0=y/a_nd;
	z0=z/a_nd;
	ax0=sma_ax/a_nd;
	ay0=sma_ay/a_nd;
	az0=sma_az/a_nd;

	
	//solving cubic
	p = (ax0*ax0 + ay0*ay0 + az0*az0) - (x0*x0 + y0*y0 + z0*z0);
	q = ( (ax0*ax0)*(ay0*ay0) + (ay0*ay0)*(az0*az0) + (az0*az0)*(ax0*ax0) - (x0*x0)*(ay0*ay0 + az0*az0) - (y0*y0)*(ax0*ax0 + az0*az0) - (z0*z0)*(ax0*ax0 + ay0*ay0) );
	r = (ax0*ax0)*(ay0*ay0)*(az0*az0) - (x0*x0)*(ay0*ay0)*(az0*az0) - (ax0*ax0)*(y0*y0)*(az0*az0) - (ax0*ax0)*(ay0*ay0)*(z0*z0);
	
	for(i=1;i<3;i++) sol_cubic[i]=0.;
	
	gsl_poly_solve_cubic(p,q,r, &sol_cubic[0],&sol_cubic[1],&sol_cubic[2]);
	
	for(i=0;i<=2;i++) {
		if (lambda < sol_cubic[i]) lambda = sol_cubic[i];
	}

	//solving semi-infinite integral
	for(i=0; i<3; i++){
		gsl_integration_workspace * workspace0 = gsl_integration_workspace_alloc (limit0);
		gsl_function F;
		
		pow_xyz[0]=0.;	pow_xyz[1]=0.;	pow_xyz[2]=0.;
		pow_xyz[i] = 1.0;
		struct my_f_params params = {ax0, ay0, az0, pow_xyz[0], pow_xyz[1], pow_xyz[2]};
		F.function = &my_f;
		F.params = &params;
		
		gsl_integration_qagiu (&F, lambda, epsabs0, epsrel0, limit0, workspace0, &result[i], &error[i]);
		
		//~ printf ("result B[%d] = % .18e\n", i+1, result[i]);
		//~ printf ("estimated error = % .18e\n", error[i]);
		//~ printf ("intervals       = %zu\n", workspace0->size);
		gsl_integration_workspace_free (workspace0);
	}
	
	Bx=result[0];By=result[1];Bz=result[2];

	return 0;

}





//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//output
/*
void heartbeat(struct reb_simulation* r){
	FILE *fpos = fopen("pos.txt","a+");
	FILE *fvel = fopen("vel.txt","a+");
	FILE *facc = fopen("acc.txt","a+");
	FILE *fprop = fopen("prop.txt","a+");

	double omega_cb=(2.*M_PI/rot_period);
    const int N = r->N;
    double PE_all=0., KE_all=0., C_all=0.;
    double Lx=0;
    double Ly=0;
    double Lz=0;
    double L_all=0;
    
	for (int i=0;i<N;i++){
		double x = r->particles[i].x; double y = r->particles[i].y; double z = r->particles[i].z;
		double vx = r->particles[i].vx; double vy = r->particles[i].vy; double vz = r->particles[i].vz;
		double ax = r->particles[i].ax; double ay = r->particles[i].ay; double az = r->particles[i].az;
		double r_ptcl = sqrt(x*x + y*y + z*z);double v_ptcl = sqrt(vx*vx + vy*vy + vz*vz);double a_ptcl = sqrt(ax*ax + ay*ay + az*az);
		//~ double cb_mass=(4./3.)*M_PI*sma_ax*sma_ay*sma_az*rho_cb;
		
		int_ellipsoid_V(x, y, z);
		PE_all = PE_all + V_ell;//total potential energy
		KE_all = KE_all + (v_ptcl*v_ptcl/2.);//kinetic energy in rotating reference frame
		C_all = C_all + ( (v_ptcl*v_ptcl/2.) + V_ell - (omega_cb*omega_cb*(x*x + y*y)/2.) );//Jacobi constant/constamt of motion in rotating reference frame
		Lx = Lx + y*vz - z*vy; Ly = Ly + z*vx - x*vz;	Lz = Lz + x*vy - y*vx;		
		
		fprintf(fpos, "%d\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",i+1,r->t, r_ptcl, x, y, z);
		fprintf(fvel, "%d\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",i+1,r->t, v_ptcl, vx, vy, vz);
		fprintf(facc, "%d\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",i+1,r->t, a_ptcl, ax, ay, az);
		
   }
   
	L_all = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
    fprintf(fpos, "\n");
    fprintf(fvel, "\n");
    fprintf(facc, "\n");
    fprintf(fprop, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",r->t, L_all, PE_all, KE_all, C_all );

    printf("Current time: %f (%f %%) \n", r->t, (100.*r->t/tmax));

	fclose(fpos);
    fclose(fvel);
    fclose(facc);
    fclose(fprop);
}
*/


void heartbeat(struct reb_simulation* r){
	FILE *fpos = fopen("pos.txt","a+");
	FILE *fvel = fopen("vel.txt","a+");
	FILE *facc = fopen("acc.txt","a+");
	FILE *fprop = fopen("prop.txt","a+");

	double omega_cb=(2.*M_PI/rot_period);
    const int N = r->N;
    double PE_all=0., KE_all=0., C_all=0.;
    double Lx=0;
    double Ly=0;
    double Lz=0;
    double L_all=0;
    
	for (int i=0;i<N;i++){
		double x = r->particles[i].x; double y = r->particles[i].y; double z = r->particles[i].z;
		double vx = r->particles[i].vx; double vy = r->particles[i].vy; double vz = r->particles[i].vz;
		double ax = r->particles[i].ax; double ay = r->particles[i].ay; double az = r->particles[i].az;
		double r_ptcl = sqrt(x*x + y*y + z*z);double v_ptcl = sqrt(vx*vx + vy*vy + vz*vz);double a_ptcl = sqrt(ax*ax + ay*ay + az*az);
		//~ double cb_mass=(4./3.)*M_PI*sma_ax*sma_ay*sma_az*rho_cb;
		
		V_ell=0;//int_ellipsoid_V(x, y, z);
		PE_all = PE_all + V_ell;//total potential energy
		KE_all = KE_all + (v_ptcl*v_ptcl/2.);//kinetic energy in rotating reference frame
		C_all = C_all + ( (v_ptcl*v_ptcl/2.) + V_ell - (omega_cb*omega_cb*(x*x + y*y)/2.) );//Jacobi constant/constamt of motion in rotating reference frame
		Lx = Lx + y*vz - z*vy; Ly = Ly + z*vx - x*vz;	Lz = Lz + x*vy - y*vx;		
		
		fprintf(fpos, "%d\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",i+1,r->t, r_ptcl, x, y, z);
		fprintf(fvel, "%d\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",i+1,r->t, v_ptcl, vx, vy, vz);
		fprintf(facc, "%d\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",i+1,r->t, a_ptcl, ax, ay, az);
		
   }
   
	L_all = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
    fprintf(fpos, "\n");
    fprintf(fvel, "\n");
    fprintf(facc, "\n");
    fprintf(fprop, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",r->t, L_all, PE_all, KE_all, C_all );

    printf("Current time: %f (%f %%) \n", r->t, (100.*r->t/tmax));

	fclose(fpos);
    fclose(fvel);
    fclose(facc);
    fclose(fprop);
}






//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//function form for gravitational potential to triaxial ellipsoid at an external point
struct my_f_params2 { double AX0; double AY0; double AZ0; double X0; double Y0; double Z0;};

double my_f2 (double u, void * p) {
	
	struct my_f_params2 * params2 = (struct my_f_params2 *)p;
	double AX0 = (params2->AX0);
	double AY0 = (params2->AY0);
	double AZ0 = (params2->AZ0);
	double X0 = (params2->X0);
	double Y0 = (params2->Y0);
	double Z0 = (params2->Z0);

	return ( AX0*AY0*AZ0*( (X0*X0/(AX0*AX0 + u)) +(Y0*Y0/(AY0*AY0 + u)) +(Z0*Z0/(AZ0*AZ0 + u)) - 1.)/( pow(((AX0*AX0 + u)*(AY0*AY0 + u)*(AZ0*AZ0 + u)), 0.5) ) );

}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//To calculate potential due to central body
int int_ellipsoid_V (double x, double y, double z) {
	
	//non-dimesionalized parameters: to avoid numerical issues with operations on large numbers
	double x0, y0, z0;
	double ax0, ay0, az0;
	
	double a_nd;

	double p, q, r;  //cubic eqn's coeff.
	double	sol_cubic[3];
	double	lambda=0.;
	double result;
	double error;
	size_t limit0=10000;
	int i;
	
	a_nd=sma_ax;
	
	x0=x/a_nd;
	y0=y/a_nd;
	z0=z/a_nd;
	ax0=sma_ax/a_nd;
	ay0=sma_ay/a_nd;
	az0=sma_az/a_nd;

	
	//solving cubic
	p = (ax0*ax0 + ay0*ay0 + az0*az0) - (x0*x0 + y0*y0 + z0*z0);
	q = ( (ax0*ax0)*(ay0*ay0) + (ay0*ay0)*(az0*az0) + (az0*az0)*(ax0*ax0) - (x0*x0)*(ay0*ay0 + az0*az0) - (y0*y0)*(ax0*ax0 + az0*az0) - (z0*z0)*(ax0*ax0 + ay0*ay0) );
	r = (ax0*ax0)*(ay0*ay0)*(az0*az0) - (x0*x0)*(ay0*ay0)*(az0*az0) - (ax0*ax0)*(y0*y0)*(az0*az0) - (ax0*ax0)*(ay0*ay0)*(z0*z0);
	
	for(i=1;i<3;i++) sol_cubic[i]=0.;
	
	gsl_poly_solve_cubic(p,q,r, &sol_cubic[0],&sol_cubic[1],&sol_cubic[2]);
	
	for(i=0;i<=2;i++) {
		if (lambda < sol_cubic[i]) lambda = sol_cubic[i];
	}

	//solving semi-infinite integral
		gsl_integration_workspace * workspace0 = gsl_integration_workspace_alloc (limit0);
		gsl_function F2;
		
		struct my_f_params2 params2 = {ax0, ay0, az0, x0, y0, z0};
		F2.function = &my_f2;
		F2.params = &params2;
		
		gsl_integration_qagiu (&F2, lambda, epsabs0, epsrel0, limit0, workspace0, &result, &error);
		
		V_ell = M_PI*G*rho_cb*a_nd*a_nd*result;
		
		gsl_integration_workspace_free (workspace0);

	return 0;

}
