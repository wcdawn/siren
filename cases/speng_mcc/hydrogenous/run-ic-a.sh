#!/bin/bash
#
#SBATCH -J PRISM_IC 
#SBATCH -e error.dat
#SBATCH -o info.dat
#SBATCH -p gradq
#SBATCH -N 2
#SBATCH --mem 4000

#cd $SLURM_SUBMIT_DIR

cwd=`pwd`
# ------------------------------------------------------------------------------ libraries
mccdir=/home/wcdawn/projects/mcc
libdir=$mccdir/lib
lib_mcclibdir=$libdir/lib.mcc.e70
lib_pwlibdir=$libdir/lib.pw.200.e70
lib_gammalibdir=$libdir/lib.gamma.e70
exe_mcc3=$mccdir/src/mcc3.x
# ------------------------------------------------------------------------------ cases

CASE=abr-ic-a

echo "Running case $CASE"

echo "executable:"
ls -Fl $exe_mcc3

rm -f input input.tmp

cp $CASE.mcc input

TFUEL='1200.0'
TCLAD='1200.0'
TCOOL='1200.0'
sed -i "s/TFUEL/${TFUEL}/" input
sed -i "s/TCLAD/${TCLAD}/" input
sed -i "s/TCOOL/${TCOOL}/" input

$exe_mcc3
mv output               $CASE.out
echo "finished"
