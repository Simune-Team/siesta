#!/bin/sh
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -v PATH
set -x
$1
# Usage:
#           qsub /path/to/sge_run.sh  command

