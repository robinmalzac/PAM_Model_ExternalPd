#include "m_pd.h"
#include <math.h>
#include "stdlib.h"
#include <stdio.h> 
#include <complex.h>
#include <string.h>

#define PI 3.14159265358979323846

#define DISCR 1000 // size of arrays in impedance_cyl
#define NMODES 3 // number of modes

static t_class *modelEul_tilde_class;

typedef struct // structure to hold the returns of impedance_cyl
{
    double Fn [DISCR];
    double Yn [DISCR];
    double fn [DISCR];
    //double complex Ze [DISCR];
}cylinder_params;

cylinder_params impedance_cyl(double, double);

typedef struct _modelEul_tilde 
{
	t_object x_obj;
	t_float gamma;
	t_float zeta;
	t_float gamma_prev; // for interpolation
	t_float zeta_prev;
	t_int fe; // sampling frequency
	t_float L; // length of resonator
	t_float R; // radius of resonator
	double p1 [NMODES]; // pn(t)
	double p2 [NMODES]; // pn'(t)
	double p1tot; // p(t) total
	double p2tot; // p'(t) total
	cylinder_params Cyl; // returns from impedance_cyl
	t_inlet *x_in2;
	t_outlet *x_out;
} t_modelEul_tilde;

void modelEul_tilde_set_gamma(t_modelEul_tilde *x, t_floatarg f)
{
	x->gamma = f;
}

t_int *modelEul_tilde_perform(t_int *w)
{
	t_modelEul_tilde *x = (t_modelEul_tilde *)(w[1]);
	t_sample *out = (t_sample *)(w[2]);
	int n = (int)(w[3]);

    double Te = 1./x->fe;
    double p1_old [NMODES];
    for(int k=0 ; k<NMODES ; k++) {
    	p1_old[k] = x->p1[k];
    }
    t_float gamma_interp, zeta_interp;

	for (int i=0; i<n; i++){
		gamma_interp = x->gamma_prev + (double)i/n*(x->gamma - x->gamma_prev); // interpolation
		zeta_interp = x->zeta_prev + (double)i/n*(x->zeta - x->zeta_prev); // interpolation

		// calculation of parameters with current gamma, zeta
    	double A = zeta_interp*(3*gamma_interp-1)/(2*sqrt(gamma_interp));
    	double B = -zeta_interp*(3*gamma_interp+1)/(8*pow(gamma_interp,3./2));
    	double C = -zeta_interp*(gamma_interp+1)/(16*pow(gamma_interp,5./2));

    	x->p1tot = 0;
    	x->p2tot = 0;

    	for(int j=0 ; j<NMODES ; j++) {
    		x->p1tot += x->p1[j];
    		x->p2tot += x->p2[j];
    	}

    	for(int j=0 ; j<NMODES ; j++) {
	        x->p1[j] = x->p1[j] + Te*x->p2[j];
	        x->p2[j] = x->p2[j] - Te*(x->Cyl.Fn[j]*(x->Cyl.Yn[j]*x->p2[j]
	        										-x->p2tot*(A+2*B*x->p1tot+3*C*pow(x->p1tot,2)))
	        										+pow(2*PI*x->Cyl.fn[j],2)*p1_old[j]);
	        p1_old[j] = x->p1[j];
	    }

        out[i] = x->p1tot/NMODES;
	}

	x->gamma_prev = x->gamma;
	x->zeta_prev = x->zeta;

	return (w+4);
}

void modelEul_tilde_dsp(t_modelEul_tilde *x, t_signal **sp)
{
	dsp_add(modelEul_tilde_perform, 3, x,
	sp[0]->s_vec, sp[0]->s_n);
}

void modelEul_tilde_free(t_modelEul_tilde *x)
{
	inlet_free(x->x_in2);
	outlet_free(x->x_out);
}

void *modelEul_tilde_new(t_floatarg f1, t_floatarg f2)
{
	t_modelEul_tilde *x = (t_modelEul_tilde *)pd_new(modelEul_tilde_class);
	x->gamma = f1;
	x->zeta = f2;
	x->gamma_prev = x->gamma;
	x->zeta_prev = x->zeta;
	x->fe = 44100;
	x->L = 0.660;
	x->R = 0.007;
	x->Cyl = impedance_cyl(x->L, x->R); // calculates resonator parameters

    char str[80];
    sprintf(str, "%f", x->Cyl.fn[0]); // print first eigenmode
	post(str);

	// initialization
	for(int i=0 ; i < NMODES ; ++i) {
		x->p1[i] = 0.0001;
		x->p2[i] = 0;
	}
	x->p1tot = 0;
	x->p2tot = 0;

	x->x_in2 = floatinlet_new(&x->x_obj, &x->zeta);
	x->x_out = outlet_new(&x->x_obj, &s_signal);
	return (void *)x;
}

void modelEul_tilde_setup(void)
{
	modelEul_tilde_class = class_new(gensym("modelEul~"),
	(t_newmethod)modelEul_tilde_new,
	0, sizeof(t_modelEul_tilde),
	CLASS_DEFAULT,
	A_DEFFLOAT,
	A_DEFFLOAT, 0);

	class_addfloat(modelEul_tilde_class,
	(t_method)modelEul_tilde_set_gamma);

	class_addmethod(modelEul_tilde_class,
	(t_method)modelEul_tilde_dsp, gensym("dsp"), 0);
}

/* Calculations below */

cylinder_params impedance_cyl(double L, double R)
{
	double c = 340;
    double rho = 1.125;
    double S = PI*pow(R,2);
    double Zc = rho*c/S; // impedance caracteristique
    int size_w = DISCR; // size of most arrays
    int start = 2*PI*10;
    int end = 2*PI*3000;
    double w[size_w], k[size_w], Zabs[size_w];
    double complex Zs[size_w], GAMMA[size_w], Ze[size_w];

    for(int i=0 ; i < size_w ; ++i) {
    	w[i] = start + i*(end-start)/(size_w-1);
    	k[i] = w[i]/c;
    	Zs[i] = Zc*(I*k[i]*0.6*R + 0.25*(k[i]*pow(R,2))); // Impedance de rayonnement
    	GAMMA[i] = I*w[i]/c + (1+I)*3e-5*sqrt(w[i]/(2*PI))/R; // amortissement et dissipation
    	Ze[i] = Zc*ctanh(GAMMA[i]*L + catanh(Zs[i]/Zc)); // Impedance d'entree
    	Zabs[i] = cabs(Ze[i]);
    }


    /* Find peaks */
	int peaks[DISCR] = { 0 };
	int pe = 0; // number of peaks

	for (int i = 1; i < size_w-1; i++) {
		if (Zabs[i]>Zabs[i-1] && Zabs[i]>Zabs[i+1])
	    {
	    	peaks[pe++] = i; // stores all peak indices
	    }
	}

	int ind[pe];
    
    for(int i=0 ; peaks[i]!=0 ; ++i) {
    	ind[i] = peaks[i];
    }

    double ZMn[pe], fn[pe], Fn[pe];
    double Yn[pe];

    for(int i=0 ; i < pe ; ++i) {
    	ZMn[i] = Zabs[ind[i]];
    	Yn[i] = 1./ZMn[i];
    	fn[i] = w[ind[i]]/(2*PI);
    	Fn[i] = 2*c/L;
    }

    // Pass the results to a structure
    cylinder_params result;
    memcpy(result.Fn, Fn, pe*sizeof(double));
    memcpy(result.Yn, Yn, pe*sizeof(double));
    memcpy(result.fn, fn, pe*sizeof(double));
    //memcpy(result.Ze, Ze, DISCR*sizeof(double complex));
    return result;
}

