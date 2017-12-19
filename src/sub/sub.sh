#MSUB -N setupfiles
#MSUB -q thompson
#MSUB -j oe
#MSUB -d ./
#MSUB -l nodes=1:ppn=10:ib,mem=5gb,walltime=24:00:00

module load legacy
module load intel_compiler_2016
module load intel_mpi_intel64/5.1.2.150
module load lammps/20151207

echo Time is `date`
echo Directory is `pwd`

bash setup_files

echo Ending Time is `date`
exit 0
