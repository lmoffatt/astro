#!/bin/bash


extract_par_beta_1(){
cd $1
a=($ls *par*.000*)
cp $a par_header

grep $model_[0-9]*[[:space:]][0-9]*[[:space:]][-.0-9]*[[:space:]][0-9]**[[:space:]][0-9]*[[:space:]][-.0-9]*[[:space:]][.0-9]**[[:space:]][0-9]*[[:space:]][-.0-9]*  *par.* -h >total

cat par_header total>../beta_par_$1.txt
rm total
cd ..
}
extract_all(){
a=($ls ./m*/*$1.*.000*)
cp $a $1_header
cat $1_header ./m*/*$1.*.00[1-9]* ./m*/*$1.*.0[1-9][0-9]* >total_$1.txt
}


extract_eq(){
a=($ls ./m*/*$1*.000*)

cp $a $1_header
cat $1_header ./m*/*$1.*.01[7-9]* ./m*/*$1.*.0[2-9][0-9]* >total_eq_$1.txt
}

extract_beta_1(){
a=($ls ./m*/*$1*.000*)

cp $a $1_header

grep $model_[0-9]*[[:space:]][0-9]*[[:space:]][-.0-9]*[[:space:]][0-9]**[[:space:]][0-9]*[[:space:]][-.0-9]*[[:space:]][1]**[[:space:]][0-9]*[[:space:]][-.0-9]*  ./m*/*$1.* -h >total

cat $1_header total>beta_$1_total.txt
rm total
}

extract_all logL
extract_beta_1 fit
extract_beta_1 sim
extract_par_beta_1 m01
extract_par_beta_1 m02
extract_par_beta_1 m03
extract_par_beta_1 m10
extract_par_beta_1 m101
extract_par_beta_1 m11
extract_par_beta_1 m12
extract_par_beta_1 m13




