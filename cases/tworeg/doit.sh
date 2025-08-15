#!/bin/sh

EXEC=../../src/siren.x

BASE=tworeg_nonuniform

cp ${BASE}.inp ${BASE}_run.inp

sed -i '' 's/refine.*$/refine 0/' ${BASE}_run.inp
${EXEC} ${BASE}_run.inp | grep 'linferr ='

sed -i '' 's/refine.*$/refine 1/' ${BASE}_run.inp
${EXEC} ${BASE}_run.inp | grep 'linferr ='

sed -i '' 's/refine.*$/refine 2/' ${BASE}_run.inp
${EXEC} ${BASE}_run.inp | grep 'linferr ='

sed -i '' 's/refine.*$/refine 3/' ${BASE}_run.inp
${EXEC} ${BASE}_run.inp | grep 'linferr ='

sed -i '' 's/refine.*$/refine 4/' ${BASE}_run.inp
${EXEC} ${BASE}_run.inp | grep 'linferr ='

sed -i '' 's/refine.*$/refine 5/' ${BASE}_run.inp
${EXEC} ${BASE}_run.inp | grep 'linferr ='

sed -i '' 's/refine.*$/refine 6/' ${BASE}_run.inp
${EXEC} ${BASE}_run.inp | grep 'linferr ='

sed -i '' 's/refine.*$/refine 7/' ${BASE}_run.inp
${EXEC} ${BASE}_run.inp | grep 'linferr ='

sed -i '' 's/refine.*$/refine 8/' ${BASE}_run.inp
${EXEC} ${BASE}_run.inp | grep 'linferr ='

sed -i '' 's/refine.*$/refine 9/' ${BASE}_run.inp
${EXEC} ${BASE}_run.inp | grep 'linferr ='

sed -i '' 's/refine.*$/refine 10/' ${BASE}_run.inp
${EXEC} ${BASE}_run.inp | grep 'linferr ='

rm -f ${BASE}_run.inp
