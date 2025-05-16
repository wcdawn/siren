#!/bin/sh

EXEC=../../src/siren.x

cp tworeg.inp tworeg_run.inp

sed -i '' 's/refine.*$/refine 0/' tworeg_run.inp
${EXEC} tworeg_run.inp | grep 'keff ='

sed -i '' 's/refine.*$/refine 1/' tworeg_run.inp
${EXEC} tworeg_run.inp | grep 'keff ='

sed -i '' 's/refine.*$/refine 2/' tworeg_run.inp
${EXEC} tworeg_run.inp | grep 'keff ='

sed -i '' 's/refine.*$/refine 3/' tworeg_run.inp
${EXEC} tworeg_run.inp | grep 'keff ='

sed -i '' 's/refine.*$/refine 4/' tworeg_run.inp
${EXEC} tworeg_run.inp | grep 'keff ='

sed -i '' 's/refine.*$/refine 5/' tworeg_run.inp
${EXEC} tworeg_run.inp | grep 'keff ='

sed -i '' 's/refine.*$/refine 6/' tworeg_run.inp
${EXEC} tworeg_run.inp | grep 'keff ='

sed -i '' 's/refine.*$/refine 7/' tworeg_run.inp
${EXEC} tworeg_run.inp | grep 'keff ='

sed -i '' 's/refine.*$/refine 8/' tworeg_run.inp
${EXEC} tworeg_run.inp | grep 'keff ='

sed -i '' 's/refine.*$/refine 9/' tworeg_run.inp
${EXEC} tworeg_run.inp | grep 'keff ='

sed -i '' 's/refine.*$/refine 10/' tworeg_run.inp
${EXEC} tworeg_run.inp | grep 'keff ='

rm -f tworeg_run.inp
