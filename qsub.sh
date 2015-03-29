#! /bin/bash
echo cd `pwd` \; "$@" | qsub -h -l mem=4gb,walltime=48:00:00
