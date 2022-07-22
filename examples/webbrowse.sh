#! /bin/sh 

# make sure it's run from the right place

if [ $# -lt 2 ] 
then 
 echo usage: examples/webbrowse.sh station run \[datadir=/data/rootified\] [port=12345]
 exit 1
fi 


if [ "$(basename $(pwd))" != "mattak" ] 
then 
  echo "Run this script from mattak (i.e. as examples/webbrowse.sh" 
  echo "Or, even better, fix it so it can run from other places :" 
  exit 1
fi 

STATION=$1 
RUN=$2 
DATA=${3-/data/rootified}
PORT=${4-12345} 

LD_LIBRARY_PATH+=:/usr/local/lib root -b examples/webbrowse.C\($STATION,$RUN,\"$DATA\",$PORT\); 

