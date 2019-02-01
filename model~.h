#ifndef _MODEL_
# define _MODEL_

#include "m_pd.h"


static t_class *model_tilde_class;
typedef struct _model_tilde {
  t_object x_obj;

  t_sample f;
  t_outlet *x_out; // signal output
  t_inlet *in2; // input of the second signal

  float ti;
  float freq;

  float shapeWidth;
  float bypass;

} t_model_tilde;

void *model_tilde_new(t_symbol *s, int argc, t_atom *argv);
void model_tilde_free(t_model_tilde *x);
void model_tilde_setup(void);

t_int *model_tilde_perform(t_int *w);
void model_tilde_dsp(t_model_tilde *x, t_signal **sp);

#endif
