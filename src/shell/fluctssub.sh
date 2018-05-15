#MSUB -N flucts_array
#MSUB -q laird
#MSUB -j oe
#MSUB -d ./
#MSUB -l nodes=1:ppn=2:intel,mem=30gb,walltime=30:00:00
#MSUB -t 0-6


cd $PBS_O_WORKDIR

echo Time is `date` > array_$MOAB_JOBARRAYINDEX.o
echo Directory is `pwd` >> array_$MOAB_JOBARRAYINDEX.o
