#! /bin/bash
root_path=../cv-c-code
num_factors=200
#root_path=../data/arxiv/cv
#num_factors=100

for i in `seq 1 5`
do
#  ./qsub.sh ./ctr --directory $root_path/cv-cf-$i --user $root_path/cf-train-$i-users.dat --item \
#  $root_path/cf-train-$i-items.dat --a 1 --b 0.01 --lambda_u 0.01 --lambda_v 0.01 \
#  --random_seed 33333 --num_factors $num_factors --save_lag 20

  for type in ofm cf
  do
  ./qsub.sh ./ctr --directory $root_path/cv-ctr-$i-$type --user $root_path/$type-train-$i-users.dat --item \
  $root_path/$type-train-$i-items.dat --a 1 --b 0.01 --lambda_u 0.01 --lambda_v 100 \
  --mult $root_path/mult.dat --theta_init $root_path/theta-vector.dat \
  --beta_init $root_path/final.beta --num_factors $num_factors --save_lag 20 --theta_opt
  done

done
