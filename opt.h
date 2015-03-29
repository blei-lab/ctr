#ifndef OPT_H
#define OPT_H
#include <gsl/gsl_vector.h>

void optimize_simplex(const gsl_vector* gamma, const gsl_vector* v, double lambda, gsl_vector* opt_x);
double f_simplex(const gsl_vector* gamma, const gsl_vector* v, double lambda, const gsl_vector* opt_x);
void df_simplex(const gsl_vector* gamma, const gsl_vector* v, double lambda, const gsl_vector* opt_x, gsl_vector* g);
bool is_feasible(const gsl_vector* x);
void simplex_projection(const gsl_vector* x, gsl_vector* x_proj, double z=1.0);

#endif // OPT_H

