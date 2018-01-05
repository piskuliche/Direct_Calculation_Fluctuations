#MSUB -N msd_calc
#MSUB -q sixhour
#MSUB -j oe
#MSUB -d ./
#MSUB -l nodes=1:ppn=2:ib,mem=5gb,walltime=6:00:00
#MSUB -t AAA-BBB

SEP=1000
START=1001000
CUR=$(( MOAB_JOBARRAYINDEX*SEP - SEP + START ))
cd $PBS_O_WORKDIR
cd FILES/$CUR

echo Time is `date` > array_$MOAB_JOBARRAYINDEX.o
echo Directory is `pwd` >> array_$MOAB_JOBARRAYINDEX.o

cp ../../msd_rot_calc ./
./msd_rot_calc > tmp.msd
 
echo Ending Time is `date` >> array_$MOAB_JOBARRAYINDEX.o
exit 0
