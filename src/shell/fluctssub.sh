#MSUB -N FLUCTSCALC
#MSUB -q sixhour
#MSUB -d ./
#MSUB -j oe
#MSUB -l nodes=1:ppn=10:intel,mem=5gb,walltime=6:00:00


module load Python/2.7
