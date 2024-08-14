#!/bin/sh

for ((INDEX = 0; INDEX < 20; INDEX++))
do

mkdir dir${INDEX}

cp listareq28_${INDEX}.txt dir${INDEX} 

cd dir${INDEX}

split -l 556 -d listareq28_${INDEX}.txt

mkdir xa0 

mkdir xa1

cp x00 xa0/x00
cp x01 xa1/x01

cd ..

cp runanauno dir${INDEX}/xa0/runanauno
cp NucleiAna.json dir${INDEX}/xa0/NucleiAna.json
cp pCentr.txt dir${INDEX}/xa0/pCentr.txt

cp runanauno dir${INDEX}/xa1/runanauno
cp NucleiAna.json dir${INDEX}/xa1/NucleiAna.json
cp pCentr.txt dir${INDEX}/xa1/pCentr.txt

sed -e "s/listTEST/x00/" RunAnalyses.C > dir${INDEX}/xa0/RunAnalyses.C

sed -e "s/listTEST/x01/" RunAnalyses.C > dir${INDEX}/xa1/RunAnalyses.C
 
cd dir${INDEX}/xa0/

sbatch runanauno

cd ../xa1

sbatch runanauno

cd ../..

done
