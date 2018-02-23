#!/usr/bin/env zsh
 
### Job name
#BSUB -J JOBNAME
 
### File / path where STDOUT & STDERR will be written
###    %J is the job ID, %I is the array ID
#BSUB -o job.%J.%I
 
### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 80 minutes you could also use this: 1:20
#BSUB -W 1:00
 
### Request memory you need for your job in TOTAL in MB
#BSUB -M 2048

### Request the number of compute slots you want to use
#BSUB -n 1
 
### Use esub for OpenMP/shared memeory jobs
#BSUB -a openmp

### Request special machine type (Intel Westmere X5675)
#BSUB -m mpi-s

### Request the node exclusivly
#BSUB -x

### Change to the work directory
cd /home/user/workdirectory
 
### Execute your application
./output 376 561 477 532 425 590
