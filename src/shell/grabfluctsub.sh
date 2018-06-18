#MSUB -N grabflucts_array
#MSUB -q sixhour
#MSUB -j oe
#MSUB -d ./
#MSUB -l nodes=1:ppn=2:intel,mem=30gb,walltime=6:00:00
#MSUB -t 0-3


cd $PBS_O_WORKDIR

echo Time is `date` > array_$MOAB_JOBARRAYINDEX.o
echo Directory is `pwd` >> array_$MOAB_JOBARRAYINDEX.o

python grab_flucts.py flucts.inp $MOAB_JOBARRAYINDEX
