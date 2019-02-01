#include "m_pd.h"
#include "math.h"

#define PI 3.14159265358979323846

static t_class *harmony_tilde_class;

typedef struct _harmony_tilde {
	t_object x_obj;
	t_sample freq;
	t_sample harm;
	t_sample ti;
	t_inlet *x_in2;
	t_outlet *x_out;
} t_harmony_tilde;

void harmony_tilde_set_freq(t_harmony_tilde *x, t_floatarg f)
{
	x->freq = f;
}

t_int *harmony_tilde_perform(t_int *w)
{
	t_harmony_tilde *x = (t_harmony_tilde *)(w[1]);
	t_sample *out = (t_sample *)(w[2]);
	int n = (int)(w[3]);

	for (int i = 0; i < n; i++){
		out[i] = 0.25*(sin(2*PI * x->freq * x->ti/44000) +
				 x->harm * (sin(4*PI * x->freq * x->ti/44000) +
				 			sin(6*PI * x->freq * x->ti/44000) +
				 			sin(8*PI * x->freq * x->ti/44000)));
		x->ti += 1;
	}

	return (w+4);
}

void harmony_tilde_dsp(t_harmony_tilde *x, t_signal **sp)
{
	dsp_add(harmony_tilde_perform, 3, x,
	sp[0]->s_vec, sp[0]->s_n);
}

void harmony_tilde_free(t_harmony_tilde *x)
{
	inlet_free(x->x_in2);
	outlet_free(x->x_out);
}

void *harmony_tilde_new(t_floatarg f1, t_floatarg f2)
{
	t_harmony_tilde *x = (t_harmony_tilde *)pd_new(harmony_tilde_class);
	x->freq = f1;
	x->harm = f2;
	x->ti = 0;
	x->x_in2 = floatinlet_new(&x->x_obj, &x->harm);
	x->x_out = outlet_new(&x->x_obj, &s_signal);
	return (void *)x;
}

void harmony_tilde_setup(void) {
	harmony_tilde_class = class_new(gensym("harmony~"),
	(t_newmethod)harmony_tilde_new,
	0, sizeof(t_harmony_tilde),
	CLASS_DEFAULT,
	A_DEFFLOAT,
	A_DEFFLOAT, 0);

	class_addfloat(harmony_tilde_class,
	(t_method)harmony_tilde_set_freq);

	class_addmethod(harmony_tilde_class,
	(t_method)harmony_tilde_dsp, gensym("dsp"), 0);
}

