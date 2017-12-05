#MSUB -N msd_rot_calc
#MSUB -q thompson
#MSUB -j oe
#MSUB -M piskuliche@ku.edu
#MSUB -m ae
#MSUB -d ./
#MSUB -l nodes=1:ppn=1:ib,mem=5gb,walltime=24:00:00

echo Time is `date`
echo Directory is `pwd`

./msd_rot_calc > tmp.msd

echo Ending Time is `date`
exit 0
