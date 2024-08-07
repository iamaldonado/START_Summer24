#!/bin/sh

for ((INDEX = 0; INDEX < 20; INDEX++))
do

mkdir dir${INDEX}

cp runanauno dir${INDEX}/runanauno
cp NucleiAna.json dir${INDEX}/NucleiAna.json
cp pCentr.txt dir${INDEX}/pCentr.txt
cp listareq28_${INDEX}.txt dir${INDEX} 

sed -e "s/listTEST/listareq28_${INDEX}/" RunAnalyses.C > dir${INDEX}/RunAnalyses.C

cd dir${INDEX}
sbatch runanauno
cd ..


done

