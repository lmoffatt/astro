#!/bin/bash


extract_par_beta_1(){
cd $1
a=($ls *par*.000*)
cp $a par_header

grep $model_[0-9]*[[:space:]][0-9]*[[:space:]][-.0-9]*[[:space:]][0-9]*[[:space:]][0-9]*[[:space:]][-.0-9]*[[:space:]][1][[:space:]][0-9]*[[:space:]][-.0-9]*  *par.* -h >temporaryyy

cat par_header temporaryyy>../beta_par_$1.txt
rm temporaryyy
cd ..
echo $1 parameters beta=1
echo $(ls beta_par_$1.txt -lh)
}


extract_par_eq_beta_1(){
cd $1
a=($ls *par*.000*)
cp $a par_header

grep $model_[0-9]*[[:space:]][0-9]*[[:space:]][-.0-9]*[[:space:]][0-9]*[[:space:]][5-9][0-9][0-9][[:space:]][-.0-9]*[[:space:]][1][[:space:]][0-9]*[[:space:]][-.0-9]*  *par.* -h >temporaryyy

cat par_header temporaryyy>../beta_par_eq_$1.txt
rm temporaryyy
cd ..
echo $1 parameters nsamples>500 beta=1
echo $(ls beta_par_eq_$1.txt -lh)
}


extract_all(){
a=($ls ./m*/*$1.*.000*)
cp $a $1_header
cat $1_header ./m*/*$1.*.00[1-9]* ./m*/*$1.*.0[1-9][0-9]* >total_$1.txt
echo all $1
echo $(ls total_$1.txt -lh)
}


extract_eq(){
a=($ls ./m*/*$1*.000*)

cp $a $1_header
cat $1_header ./m*/*$1.*.01[7-9]* ./m*/*$1.*.0[2-9][0-9]* >total_eq_$1.txt
echo equilibrio  $1
echo $(ls total_eq_$1.txt -lh)

}

extract_beta_eq_1(){
a=($ls ./m*/*$1*.000*)
cp $a $1_header
grep $model_[0-9]*[[:space:]][0-9]*[[:space:]][-.0-9]*[[:space:]][0-9]**[[:space:]][5-9][0-9][0-9][[:space:]][-.0-9]*[[:space:]][1][[:space:]][0-9]*[[:space:]][-.0-9]*  ./m*/*$1.* -h >temporaryyy

cat $1_header temporaryyy>beta_eq_$1_total.txt
rm temporaryyy
echo beta=1 nsamples>500 $1
echo $(ls beta_eq_$1_total.txt -lh)

}
extract_beta_1(){
a=($ls ./m*/*$1*.000*)
cp $a $1_header
grep $model_[0-9]*[[:space:]][0-9]*[[:space:]][-.0-9]*[[:space:]][0-9]*[[:space:]][0-9]*[[:space:]][-.0-9]*[[:space:]][1][[:space:]][0-9]*[[:space:]][-.0-9]*  ./m*/*$1.* -h >temporaryyy

cat $1_header temporaryyy>beta_$1_total.txt
rm temporaryyy
echo beta=1  $1
echo $(ls beta_$1_total.txt -lh)

}


extract_all logL
extract_beta_1 fit
extract_beta_eq_1 sim
extract_par_eq_beta_1 m01
extract_par_eq_beta_1 m02
extract_par_eq_beta_1 m03
extract_par_eq_beta_1 m10
extract_par_eq_beta_1 m101
extract_par_eq_beta_1 m11
extract_par_eq_beta_1 m12
extract_par_eq_beta_1 m13




