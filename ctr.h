// class for ctr
//
#ifndef CTR_H
#define CTR_H

#include "utils.h"
#include "corpus.h"
#include "data.h"

struct ctr_hyperparameter {
  double a;
  double b;
  double lambda_u;
  double lambda_v;
  double learning_rate;
  double alpha_smooth;
  int    random_seed;
  int    max_iter;
  int    save_lag;
  int    theta_opt;
  int    ctr_run;
  int    lda_regression;
  
  void set(double aa, double bb, 
           double lu, double lv, 
           double lr, double as,
           int rs, int mi, int sl,    
           int to, int cr, int lda_r) {
    a = aa; b = bb; 
    lambda_u = lu; lambda_v = lv; 
    learning_rate = lr;
    alpha_smooth = as;
    random_seed = rs; max_iter = mi;
    save_lag = sl; theta_opt = to;
    ctr_run = cr; lda_regression = lda_r;
  }

  void save(char* filename) {
    FILE * file = fopen(filename, "w");
    fprintf(file, "a = %.4f\n", a);
    fprintf(file, "b = %.4f\n", b);
    fprintf(file, "lambda_u = %.4f\n", lambda_u);
    fprintf(file, "lambda_v = %.4f\n", lambda_v);
    fprintf(file, "learning_rate = %.6f\n", learning_rate);
    fprintf(file, "alpha_smooth = %.6f\n", alpha_smooth);
    fprintf(file, "random seed = %d\n", (int)random_seed);
    fprintf(file, "max iter = %d\n", max_iter);
    fprintf(file, "save lag = %d\n", save_lag);
    fprintf(file, "theta opt = %d\n", theta_opt);
    fprintf(file, "ctr run = %d\n", ctr_run);
    fprintf(file, "lda_regression = %d\n", lda_regression);
    fclose(file);
  }
};

class c_ctr {
public:
  c_ctr();
  ~c_ctr();
  void read_init_information(const char* theta_init_path, 
                             const char* beta_init_path, 
                             const c_corpus* c, double alpha_smooth);

  void set_model_parameters(int num_factors, 
                            int num_users, 
                            int num_items);

  void learn_map_estimate(const c_data* users, const c_data* items, 
                          const c_corpus* c, const ctr_hyperparameter* param,
                          const char* directory);

  void stochastic_learn_map_estimate(const c_data* users, const c_data* items, 
                                     const c_corpus* c, const ctr_hyperparameter* param, 
                                     const char* directory);

  void init_model(int ctr_run);

  double doc_inference(const c_document* doc, const gsl_vector* theta_v, 
                       const gsl_matrix* log_beta, gsl_matrix* phi, 
                       gsl_vector* gamma, gsl_matrix* word_ss, 
                       bool update_word_ss);
public:
  gsl_matrix* m_beta;
  gsl_matrix* m_theta;

  gsl_matrix* m_U;
  gsl_matrix* m_V;

  int m_num_factors; // m_num_topics
  int m_num_items; // m_num_docs
  int m_num_users; // num of users
};

#endif // CTR_H
