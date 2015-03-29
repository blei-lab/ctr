#! /bin/bash
#root_path=../data/arxiv/cv
#num_factors=100

#for i in `seq 1 5`
#do
##  ./qsub.sh ./ctr --directory $root_path/cv-cf-$i --user $root_path/cf-train-$i-users.dat --item \
##  $root_path/cf-train-$i-items.dat --a 1 --b 0.01 --lambda_u 0.01 --lambda_v 0.01 \
##  --random_seed 33333 --num_factors $num_factors --save_lag 20
#
#  for type in ofm cf
#  do
#  ./qsub.sh ./ctr --directory $root_path/cv-ctr-$i-$type --user $root_path/$type-train-$i-users.dat --item \
#  $root_path/$type-train-$i-items.dat --a 1 --b 0.01 --lambda_u 0.01 --lambda_v 100 \
#  --mult $root_path/mult.dat --theta_init $root_path/theta-vector.dat \
#  --beta_init $root_path/final.beta --num_factors $num_factors --save_lag 20 --theta_opt
#  done
#
#done

#for i in 0 1 2 3 4
#do
#  for K in 200
#  do 
#    for lambda in 0.01 0.1 1 10 100 1000 5000
#    do
#      ./qsub.sh ./ctr --directory ../data/citeulike/data/cv-in-matrix/cf-fold-$i-K-$K-lambda-$lambda \
#      --user ../data/citeulike/data/cv-in-matrix/fold-$i-users.train \
#      --item ../data/citeulike/data/cv-in-matrix/fold-$i-items.train \
#      --lambda_u 0.01 --lambda_v $lambda --num_factors $K --save_lag 50 \
#      --max_iter 100
#
#      ./qsub.sh ./ctr --directory ../data/citeulike/data/cv-in-matrix/ctr-fold-$i-K-$K-lambda-$lambda \
#      --user ../data/citeulike/data/cv-in-matrix/fold-$i-users.train \
#      --item ../data/citeulike/data/cv-in-matrix/fold-$i-items.train \
#      --lambda_u 0.01 --lambda_v $lambda --num_factors $K \
#      --mult ../data/citeulike/data/mult.dat \
#      --theta_init ../data/citeulike/data/lda-$K/final.doc.states \
#      --beta_init ../data/citeulike/data/lda-$K/final.topics --num_factors $K \
#      --save_lag 50 --alpha_smooth 1 --max_iter 100 --theta_opt
#
#      ./qsub.sh ./ctr --directory ../data/citeulike/data/cv-out-of-matrix/ctr-fold-$i-K-$K-lambda-$lambda \
#      --user ../data/citeulike/data/cv-out-of-matrix/fold-$i-users.train \
#      --item ../data/citeulike/data/cv-out-of-matrix/fold-$i-items.train \
#      --lambda_u 0.01 --lambda_v $lambda --num_factors $K \
#      --mult ../data/citeulike/data/mult.dat \
#      --theta_init ../data/citeulike/data/lda-$K/final.doc.states \
#      --beta_init ../data/citeulike/data/lda-$K/final.topics --num_factors $K \
#      --save_lag 50 --alpha_smooth 1 --max_iter 100 --theta_opt
#
#      if [ "$lambda" == 10 ]; then
#      ./qsub.sh ./ctr --directory ../data/citeulike/data/cv-in-matrix/lda-fold-$i-K-$K-lambda-$lambda \
#      --user ../data/citeulike/data/cv-in-matrix/fold-$i-users.train \
#      --item ../data/citeulike/data/cv-in-matrix/fold-$i-items.train \
#      --lambda_u 0.01 --lambda_v $lambda --num_factors $K \
#      --mult ../data/citeulike/data/mult.dat \
#      --theta_init ../data/citeulike/data/lda-$K/final.doc.states \
#      --beta_init ../data/citeulike/data/lda-$K/final.topics --num_factors $K \
#      --save_lag 50 --alpha_smooth 1 --max_iter 100 --lda_regression
#
#      ./qsub.sh ./ctr --directory ../data/citeulike/data/cv-out-of-matrix/lda-fold-$i-K-$K-lambda-$lambda \
#      --user ../data/citeulike/data/cv-out-of-matrix/fold-$i-users.train \
#      --item ../data/citeulike/data/cv-out-of-matrix/fold-$i-items.train \
#      --lambda_u 0.01 --lambda_v $lambda --num_factors $K \
#      --mult ../data/citeulike/data/mult.dat \
#      --theta_init ../data/citeulike/data/lda-$K/final.doc.states \
#      --beta_init ../data/citeulike/data/lda-$K/final.topics --num_factors $K \
#      --save_lag 50 --alpha_smooth 1 --max_iter 100 --lda_regression
#    fi
#    done
#  done
#done

rootpath=../data/mendeley/
for i in 0 1 2 3 4
do
  for K in 500
  do 
    for lambda in 10000 20000 50000 
    do
      condor_run "./ctr-condor --directory $rootpath/condor-result/cv-in-matrix/ctr-fold-$i-K-$K-lambda-$lambda \
      --user $rootpath/cv-in-matrix/fold-$i-users.train \
      --item $rootpath/cv-in-matrix/fold-$i-items.train \
      --lambda_u 0.01 --lambda_v $lambda --num_factors $K \
      --mult $rootpath/mult.dat \
      --theta_init $rootpath/lda-$K/final.doc.states \
      --beta_init $rootpath/lda-$K/final.topics --num_factors $K --save_lag 2000 \
      --learning_rate 0.002 --random_seed 939384 --max_iter 400  --alpha_smooth 0.1" >> logs/ctr-in-$i-${lambda}.out &

      condor_run "./ctr-condor --directory $rootpath/condor-result/cv-out-of-matrix/ctr-fold-$i-K-$K-lambda-$lambda \
      --user $rootpath/cv-out-of-matrix/fold-$i-users.train \
      --item $rootpath/cv-out-of-matrix/fold-$i-items.train \
      --lambda_u 0.01 --lambda_v $lambda --num_factors $K \
      --mult $rootpath/mult.dat \
      --theta_init $rootpath/lda-$K/final.doc.states \
      --beta_init $rootpath/lda-$K/final.topics --num_factors $K --save_lag 2000 \
      --learning_rate 0.002 --random_seed 939384 --max_iter 400 --alpha_smooth 0.1" >> logs/ctr-out-$i-${lambda}.out &
    done
  done
done
