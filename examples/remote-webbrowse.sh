#! /bin/sh 

# make sure it's run from the right place

if [ $# -lt 2 ] 
then 
 echo usage: examples/remote-webbrowse.sh station run [datadir=/data/handcarry22/rootified] [port=12345]
 exit 1
fi 


if [ "$(basename $(pwd))" != "mattak" ] 
then 
  echo "Run this script from mattak (i.e. as examples/webbrowse.sh" 
  echo "Or, even better, fix it so it can run from other places :" 
  exit 1
fi 

if [ -z $RNO_G_PASSWORD ] ; 
then
  echo "the environmental variable RNO_G_PASSWORD needs to be defined!"
  exit 1 
fi 

#nothing to see here please move along
MOD_PASSWORD=`echo $RNO_G_PASSWORD | sed 's/@/AT/g' | sed 's/:/COLON/g'  | sed 's/!/BANG/g' | sed 's/\./DOT/g'`
 

STATION=$1 
RUN=$2 
DATA=${3-/data/handcarry22/rootified}
PORT=${4-12345} 

LD_LIBRARY_PATH+=:/usr/local/lib root -b examples/webbrowse.C\($STATION,$RUN,\"https://rno-g-alt:$MOD_PASSWORD@rno-g.uchicago.edu/$DATA\",$PORT\); 

