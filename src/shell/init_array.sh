#MSUB -N init_array
#MSUB -q laird
#MSUB -j oe
#MSUB -d ./
#MSUB -l nodes=1:ppn=2:intel,mem=30gb,walltime=20:00:00
#MSUB -t 0-99


cd $PBS_O_WORKDIR

python init_flucts.py $MOAB_JOBARRAYINDEX flucts.inp msd water
python init_flucts.py $MOAB_JOBARRAYINDEX flucts.inp c2 water

