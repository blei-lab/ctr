#include "opt.h"
#include "utils.h"
#include <gsl/gsl_sort_vector.h>

// projection gradient algorithm
void optimize_simplex(const gsl_vector* gamma, const gsl_vector* v, double lambda, gsl_vector* opt_x) {
  size_t size = opt_x->size;
  gsl_vector* g = gsl_vector_alloc(size);
  gsl_vector* x_bar = gsl_vector_alloc(size);
  gsl_vector* opt_x_old = gsl_vector_alloc(size);
  gsl_vector_memcpy(opt_x_old, opt_x); // save the old value
  double f_old = f_simplex(gamma, v, lambda, opt_x);
  //printf("f_old: %0.10f -> ", f_old);

  df_simplex(gamma, v, lambda, opt_x, g);
  double ab_sum = gsl_blas_dasum(g);
  if (ab_sum > 1.0) gsl_vector_scale(g, 1.0/ab_sum); // rescale the gradient

  gsl_blas_daxpy(-1, g, opt_x);
  simplex_projection(opt_x, x_bar);
  gsl_vector_sub(x_bar, opt_x_old);
  double r = 0;
  gsl_blas_ddot(g, x_bar, &r);
  r *= 0.5;

  double beta = 0.5;
  double f_new;
  double t = beta;
  int iter = 0;
  while(++iter < 100) {
    gsl_vector_memcpy(opt_x, opt_x_old);
    gsl_blas_daxpy(t, x_bar, opt_x);

    f_new = f_simplex(gamma, v, lambda, opt_x);
    if (f_new > f_old + r * t) t = t * beta;
    else break;
  }
  //printf("f_new %0.10f\n", f_new);

  if (!is_feasible(opt_x))  printf("sth is wrong, not feasible. you've got to check it ...\n");
  
  gsl_vector_free(g);
  gsl_vector_free(opt_x_old);
  gsl_vector_free(x_bar);
}

double f_simplex(const gsl_vector* gamma, const gsl_vector* v, 
                 double lambda, const gsl_vector* opt_x) {
  double f = 0.0, val;

  gsl_vector* y = gsl_vector_alloc(opt_x->size);
  gsl_vector_memcpy(y, opt_x);
  vct_log(y);

  gsl_blas_ddot(y, gamma, &f);

  gsl_vector_memcpy(y, v);
  gsl_vector_sub(y, opt_x);
  gsl_blas_ddot(y, y, &val);
  f -= 0.5 * lambda * val;
  gsl_vector_free(y);
  
  return -f;
}

void df_simplex(const gsl_vector* gamma, const gsl_vector* v, 
                double lambda, const gsl_vector* opt_x, 
                gsl_vector* g) {
  gsl_vector_memcpy(g, opt_x);
  gsl_vector_sub(g, v);
  gsl_vector_scale(g, -lambda);

  gsl_vector* y = gsl_vector_alloc(opt_x->size);
  gsl_vector_memcpy(y, gamma);
  gsl_vector_div(y, opt_x);
  gsl_vector_add(g, y);
  gsl_vector_scale(g, -1.0);
  gsl_vector_free(y);
}

bool is_feasible(const gsl_vector* x) {
 double val;
 double sum  = 0;
 for (size_t i = 0; i < x->size-1; i ++) {
   val = vget(x, i);
   if (val < 0 || val >1) return false;
   sum += val;
   if (sum > 1) return false;
 }
 return true;
}

// project x on to simplex (using // http://www.cs.berkeley.edu/~jduchi/projects/DuchiShSiCh08.pdf)
void simplex_projection(const gsl_vector* x, gsl_vector* x_proj, double z) {
  gsl_vector_memcpy(x_proj, x);
  gsl_sort_vector(x_proj);
  double cumsum = -z, u;
  int j = 0;
  int i; // this has to be int, not size_t
  for (i = (int)x->size-1; i >= 0; i --) {
    u = vget(x_proj, i);
    cumsum += u;
    if (u > cumsum/(j+1)) j++;
    else break;
  }
  double theta = cumsum/j;
  for (i = 0; i < (int)x->size; i ++) {
    u = vget(x, i)-theta;
    if (u <= 0) u = 0.0;
    vset(x_proj, i, u);
  }
  vnormalize(x_proj); // fix the normaliztion issue due to numerical errors
}
