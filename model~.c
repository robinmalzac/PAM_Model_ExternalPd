#include "model~.h"

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#define PI 3.14159265358979323846

// in order to silence unused mandatory parameters
#define UNUSED(x) (void)(x)

#define AUTONORM 1

void *model_tilde_new(t_symbol *s, int argc, t_atom *argv){

	t_model_tilde *x = (t_model_tilde *)pd_new(model_tilde_class);
	x->in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
	x->x_out = outlet_new(&x->x_obj, &s_signal);

	// Declare inlets
	inlet_new(&x->x_obj, &x->x_obj.ob_pd,
        gensym("float"), gensym("set_threshold"));
	inlet_new(&x->x_obj, &x->x_obj.ob_pd,
        gensym("float"), gensym("set_bypass"));

	x->freq = 440.0;
	x->ti = 0.0;

	x->shapeWidth = 1;
	x->bypass = 0;

	return (void *)x;

}

void model_tilde_set_threshold(t_model_tilde *x, t_floatarg f){

	x->freq = f;

}
void model_tilde_set_bypass(t_model_tilde *x, t_floatarg f){

	x->bypass = f;
}

void model_tilde_free(t_model_tilde *x){

	outlet_free(x->x_out);

}

void model_tilde_setup(void) {
	model_tilde_class = class_new(gensym("model~"),
        (t_newmethod)model_tilde_new,
        (t_method)model_tilde_free,
        sizeof(t_model_tilde),
        CLASS_DEFAULT,
        A_GIMME, 0);


 	class_addmethod(model_tilde_class,
        (t_method)model_tilde_dsp, gensym("dsp"), 0);


	class_addmethod(model_tilde_class,
        (t_method)model_tilde_set_bypass, gensym("set_bypass"),
        A_DEFFLOAT, 0);

	class_addmethod(model_tilde_class,
        (t_method)model_tilde_set_threshold, gensym("set_threshold"),
        A_DEFFLOAT, 0);

 	CLASS_MAINSIGNALIN(model_tilde_class, t_model_tilde, f); // permet d'utiliser
							// des signaux. Integre un inlet
}

t_int *model_tilde_perform(t_int *w){

	/* We get the object in order to access the buffer */
  	t_model_tilde *x = (t_model_tilde *)(w[1]);
	/* here is a pointer to the t_sample arrays that hold the 2 input signals */
	t_sample  *in1 =    (t_sample *)(w[2]);
	/* here is a pointer to the t_sample arrays that hold the 2 input signals */
	t_sample  *in2 =    (t_sample *)(w[3]);
	/* here comes the signalblock that will hold the output signal */
	t_sample  *out =    (t_sample *)(w[4]);
	/* all signalblocks are of the same length */
	int          n =           (int)(w[5]);



	for (int i = 0; i < n; i++){

		out[i] = sin(2*PI * x->freq * x->ti/44000);
		x->ti += 1;

	}

	/* return a pointer to the dataspace for the next dsp-object */
	return (w+6);

}
void model_tilde_dsp(t_model_tilde *x, t_signal **sp){

	dsp_add(model_tilde_perform, 5, x,
          sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n);


}
