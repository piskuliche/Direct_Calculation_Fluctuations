#MSUB -N flucts_array
#MSUB -q laird
#MSUB -j oe
#MSUB -d ./
#MSUB -l nodes=1:ppn=2:intel,mem=30gb,walltime=30:00:00


cd $PBS_O_WORKDIR

echo Time is `date` > init_array.out
echo Directory is `pwd` >> init_array.out
