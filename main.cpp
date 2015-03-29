#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "ctr.h"

gsl_rng * RANDOM_NUMBER = NULL;

void print_usage_and_exit() {
  // print usage information
  printf("*********************************collaborative topic models for recommendations************************\n");
  printf("Authors: Chong Wang, chongw@cs.princeton.edu, Computer Science Department, Princeton University.\n");
  printf("usage:\n");
  printf("      ctr [options]\n");
  printf("      --help:           print help information\n");

  printf("\n");
  printf("      --directory:      save directory, required\n");

  printf("\n");
  printf("      --user:           user file, required\n");
  printf("      --item:           item file, required\n");
  printf("      --a:              positive item weight, default 1\n");
  printf("      --b:              negative item weight, default 0.01 (b < a)\n");
  printf("      --lambda_u:       user vector regularizer, default 0.01\n");
  printf("      --lambda_v:       item vector regularizer, default 100\n");
  printf("      --learning_rate:  stochastic version for large datasets, default -1. Stochastic learning will be called when > 0.\n");
  printf("      --alpha_smooth:   alpha smooth, default [0.0]\n");
  printf("\n");

  printf("      --random_seed:    the random seed, default from the current time\n");
  printf("      --save_lag:       the saving lag, default 20 (-1 means no savings for intermediate results)\n");
  printf("      --max_iter:       the max number of iterations, default 200\n");
  printf("\n");

  printf("      --num_factors:    the number of factors, default 200\n");
  printf("      --mult:           mult file, in lda-c format, optional, if not provided, it's the matrix factorization\n");
  printf("      --theta_init:     topic proportions file from lda, optional (required if mult file is provided)\n");
  printf("      --beta_init:      topic distributions file from lda, optional (required if mult file is provided)\n");
  printf("      --theta_opt:      optimize theta or not, optional, default not\n");
  printf("      --lda_regression: run lda regression, default not\n");

  printf("*******************************************************************************************************\n");

  exit(0);
}

int main(int argc, char* argv[]) {
  if (argc < 2) print_usage_and_exit();

  char filename[500];
  int theta_opt = 0;
  int lda_regression = 0;

  const char* const short_options = "hd:x:i:a:b:u:v:r:s:m:k:t:e:y:z:w:";
  const struct option long_options[] = {
    {"help",          no_argument,       NULL, 'h'},
    {"directory",     required_argument, NULL, 'd'},
    {"user",          required_argument, NULL, 'x'},
    {"item",          required_argument, NULL, 'i'},
    {"a",             required_argument, NULL, 'a'},
    {"b",             required_argument, NULL, 'b'},
    {"lambda_u",      required_argument, NULL, 'u'},
    {"lambda_v",      required_argument, NULL, 'v'},
    {"random_seed",   required_argument, NULL, 'r'},
    {"save_lag",      required_argument, NULL, 's'},
    {"max_iter",      required_argument, NULL, 'm'},
    {"num_factors",   required_argument, NULL, 'k'},
    {"mult",          required_argument, NULL, 't'},
    {"theta_init",    required_argument, NULL, 'e'},
    {"beta_init",     required_argument, NULL, 'y'},
    {"learning_rate", required_argument, NULL, 'z'},
    {"alpha_smooth",  required_argument, NULL, 'w'},
    {"theta_opt",     no_argument, &theta_opt, 1},
    {"lda_regression",no_argument, &lda_regression, 1},
    {NULL, 0, NULL, 0}};

  char*  directory = NULL;

  char*  user_path = NULL;
  char*  item_path = NULL;
  double a = 1.0;
  double b = 0.01;
  double lambda_u = 0.01;
  double lambda_v = 100;
  double learning_rate = -1;
  double alpha_smooth = 0.0;

  time_t t; time(&t);
  long   random_seed = (long) t;
  int    save_lag = 20;
  int    max_iter = 200;

  int    num_factors = 200;
  char*  mult_path = NULL;
  char*  theta_init_path = NULL;
  char*  beta_init_path = NULL;

  int cc = 0; 
  while(true) {
    cc = getopt_long(argc, argv, short_options, long_options, NULL);
    switch(cc) {
      case 'h':
        print_usage_and_exit();
        break;
      case 'd':
        directory = optarg;
        break;
      case 'x':
        user_path = optarg;
        break;
      case 'i':
        item_path = optarg;
        break;
      case 'a':
        a = atof(optarg);
        break;
      case 'b':
        b = atof(optarg);
        break;
      case 'u':
        lambda_u = atof(optarg);
        break;
      case 'v':
        lambda_v = atof(optarg);
        break;
      case 'z':
        learning_rate = atof(optarg);
        break;
      case 'w':
        alpha_smooth = atof(optarg);
        break;
      case 'r':
        random_seed = atoi(optarg);
        break;
      case 's':
        save_lag = atoi(optarg);
        break;
      case 'm':
        max_iter =  atoi(optarg);
        break;    
      case 'k':
        num_factors = atoi(optarg);
        break;
      case 't':
        mult_path = optarg;
        break;
      case 'e':
        theta_init_path = optarg;
        break;
      case 'y':
        beta_init_path = optarg;
        break;
      case -1:
        break;
      case '?':
        print_usage_and_exit();
        break;
      default:
        break;
    }
    if (cc == -1)
      break;
  }

  /// print information
  printf("\n************************************************************************************************\n");
  
  if (!dir_exists(directory)) make_directory(directory);
  printf("result directory: %s\n", directory);

  if (!file_exists(user_path)) {
    printf("user file %s doesn't exist! quit ...\n", user_path);
    exit(-1);
  }
  printf("user file: %s\n", user_path);

  if (!file_exists(item_path)) {
    printf("item file %s doesn't exist! quit ...\n", item_path);
    exit(-1);
  }
  printf("item file: %s\n", item_path);

  printf("a: %.4f\n", a);
  printf("b: %.4f\n", b);
  printf("lambda_u: %.4f\n", lambda_u);
  printf("lambda_v: %.4f\n", lambda_v);
  printf("learning_rate: %.5f\n", learning_rate);
  printf("alpha_smooth: %.5f\n", alpha_smooth);
  printf("random seed: %d\n", (int)random_seed);
  printf("save lag: %d\n", save_lag);
  printf("max iter: %d\n", max_iter);
  printf("number of factors: %d\n", num_factors);

  if (mult_path != NULL) {
    if (!file_exists(item_path)) {
      printf("mult file %s doesn't exist! quit ...\n", mult_path);
      exit(-1);
    }
    printf("mult file: %s\n", mult_path);
      
    if (theta_init_path == NULL) {
      printf("topic proportions file must be provided ...\n");
      exit(-1);
    }
    if (!file_exists(theta_init_path)) {
      printf("topic proportions file %s doesn't exist! quit ...\n", theta_init_path);
      exit(-1);
    }
    printf("topic proportions file: %s\n", theta_init_path);

    if (beta_init_path == NULL) {
      printf("topic distributions file must be provided ...\n");
      exit(-1);
    }
    if (!file_exists(beta_init_path)) {
      printf("topic distributions file %s doesn't exist! quit ...\n", beta_init_path);
      exit(-1);
    }
    printf("topic distributions file: %s\n", beta_init_path);
    if (theta_opt) printf("theta optimization: True\n");
    else printf("theta optimization: false\n");
  }
  else if (theta_opt) {
    printf("theta optimization: false");
    printf("(theta_opt has no effect, back to default value: false)\n");
    theta_opt = 0;
  }

  printf("\n");

  /// save the settings
  int ctr_run = 1;
  if (mult_path == NULL) ctr_run = 0;
  ctr_hyperparameter ctr_param;
  ctr_param.set(a, b, lambda_u, lambda_v, learning_rate, alpha_smooth,
      random_seed, max_iter, save_lag, theta_opt, ctr_run, lda_regression);
  sprintf(filename, "%s/settings.txt", directory); 
  ctr_param.save(filename);
  
  /// init random numbe generator
  RANDOM_NUMBER = new_random_number_generator(random_seed);

  // read users
  printf("reading user matrix from %s ...\n", user_path);
  c_data* users = new c_data(); 
  users->read_data(user_path);
  int num_users = (int)users->m_vec_data.size();

  // read items
  printf("reading item matrix from %s ...\n", item_path);
  c_data* items = new c_data(); 
  items->read_data(item_path);
  int num_items = (int)items->m_vec_data.size();

  // create model instance
  c_ctr* ctr = new c_ctr();
  ctr->set_model_parameters(num_factors, num_users, num_items);

  c_corpus* c = NULL;
  if (mult_path != NULL) {
    // read word data
    c = new c_corpus();
    c->read_data(mult_path);
    ctr->read_init_information(theta_init_path, beta_init_path, c, alpha_smooth);
  }

  if (learning_rate <= 0) {
    ctr->learn_map_estimate(users, items, c, &ctr_param, directory);
  } else {
    ctr->stochastic_learn_map_estimate(users, items, c, &ctr_param, directory);
  }

  free_random_number_generator(RANDOM_NUMBER);
  if (c != NULL) delete c;

  delete ctr;
  delete users;
  delete items;
  return 0;
}
