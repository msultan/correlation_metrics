#!/bin/bash
#PBS -N kinase
#PBS -e /home/msultan/research/kinase/data/2013-2-1/g-chi-450micro/log.error
#PBS -o /home/msultan/research/kinase/data/2013-2-1/g-chi-450micro/log.out
#PBS -l nodes=2:ppn=24 
#PBS -l walltime=02:00:00

echo "Hello"
PBS_O_WORKDIR='/home/msultan/research/kinase/data/2013-2-1/g-chi-450micro/'
export PBS_O_WORKDIR
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
ipcontroller --ip="*" &
sleep 5

mpirun --hostfile $PBS_NODEFILE --bynode ipengine &
sleep 30

/home/msultan/software/epd_free-7.3-2-rh5-x86_64/bin/python2.7 $PBS_O_WORKDIR/g-chi-mutinf.py -d $PBS_O_WORKDIR -t 20 -i 1
ipcluster stop --profile=default

