#MSUB -N direct_calc_nve
#MSUB -q thompson
#MSUB -j oe
#MSUB -d ./
#MSUB -l nodes=1:ppn=10:ib,mem=5gb,walltime=24:00:00

module load legacy
module load intel_compiler_2016
module load intel_mpi_intel64/5.1.2.150
module load lammps/11Aug17

echo Time is `date`
echo Directory is `pwd`

mpirun lmp_mpi < in.nve -screen none
python set_msd_calcs.py -inp ${PWD##*/} -mol acn
./msd_rot_calc < msd_rot_calc.in
python set_msd_calcs.py -inp ${PWD##*/} -mol co2
./msd_rot_calc < msd_rot_calc.in

echo Ending Time is `date`
exit 0
