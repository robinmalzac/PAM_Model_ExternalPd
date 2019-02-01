#include "m_pd.h"
#include <math.h>
#include "stdlib.h"
#include <stdio.h> 
#include <complex.h>
#include <string.h>

#define PI 3.14159265358979323846

#define DISCR 1000 // size of arrays in impedance_cyl

static t_class *model_tilde_class;

typedef struct // structure to hold the returns of impedance_cyl
{
    double Fn [DISCR];
    double Yn [DISCR];
    double fn [DISCR];
    //double complex Ze [DISCR];
}cylinder_params;

cylinder_params impedance_cyl(double, double);

typedef struct _model_tilde 
{
	t_object x_obj;
	t_sample gamma;
	t_sample zeta;
	t_int fe; // sampling frequency
	t_int init; // initialization
	t_float L; // length of resonator
	t_float R; // radius of resonator
	cylinder_params Cyl; // returns from impedance_cyl
	t_inlet *x_in2;
	t_outlet *x_out;
} t_model_tilde;

void model_tilde_set_gamma(t_model_tilde *x, t_floatarg f)
{
	x->gamma = f;
}

t_int *model_tilde_perform(t_int *w)
{
	t_model_tilde *x = (t_model_tilde *)(w[1]);
	t_sample *out = (t_sample *)(w[2]);
	int n = (int)(w[3]);

	// calculation of parameters with current gamma, zeta
	t_float F0 = x->zeta*(1-x->gamma)*sqrt(x->gamma);
    t_float A = x->zeta*(3*x->gamma-1)/(2*sqrt(x->gamma));
    t_float B = -x->zeta*(3*x->gamma+1)/(8*pow(x->gamma,3./2));
    t_float C = -x->zeta*(x->gamma+1)/(16*pow(x->gamma,5./2));

	for (int i = 0; i < n; i++){
		if(x->init<2) {
			t_float p_0 = F0/(1-A);
			t_float dp_0 = 0;
			t_float p_1 = p_0 + dp_0/(x->fe);
			if(x->init==0) {
				out[i] = p_0;
			}
		    else if(x->init==1) {
    			out[i] = p_1;
		    }
		    x->init++;
		}
		else {
        	out[i] = 2*out[i-1] - out[i-2]
               		 -1./x->fe*x->Cyl.Fn[1]*(x->Cyl.Yn[1]-A)
               		 -2*B*out[i-2] 
               		 -3*C*pow(out[i-2],2)*(out[i-1]-out[i-2])
               		 -pow(2*PI*x->Cyl.fn[1]/x->fe,2)*out[i-2];
		}
	}

	return (w+4);
}

void model_tilde_dsp(t_model_tilde *x, t_signal **sp)
{
	dsp_add(model_tilde_perform, 3, x,
	sp[0]->s_vec, sp[0]->s_n);
}

void model_tilde_free(t_model_tilde *x)
{
	inlet_free(x->x_in2);
	outlet_free(x->x_out);
}

void *model_tilde_new(t_floatarg f1, t_floatarg f2)
{
	t_model_tilde *x = (t_model_tilde *)pd_new(model_tilde_class);
	x->gamma = f1;
	x->zeta = f2;
	x->fe = 44100;
	x->init = 0;
	x->L = 0.660;
	x->R = 0.015;
	x->Cyl = impedance_cyl(x->L, x->R); // calculates resonator parameters
	x->x_in2 = floatinlet_new(&x->x_obj, &x->zeta);
	x->x_out = outlet_new(&x->x_obj, &s_signal);
	return (void *)x;
}

void model_tilde_setup(void)
{
	model_tilde_class = class_new(gensym("model~"),
	(t_newmethod)model_tilde_new,
	0, sizeof(t_model_tilde),
	CLASS_DEFAULT,
	A_DEFFLOAT,
	A_DEFFLOAT, 0);

	class_addfloat(model_tilde_class,
	(t_method)model_tilde_set_gamma);

	class_addmethod(model_tilde_class,
	(t_method)model_tilde_dsp, gensym("dsp"), 0);
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

