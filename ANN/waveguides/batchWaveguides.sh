#!/usr/bin/env

echo "Bash version ${BASH_VERSION}..."

numW=1
numH=1



for iw in $(seq 0 1 $numW);
do
  for ih in $(seq 0 1 $numH);
#for iw in $W;
  do
    python3 genData.py ${iw} ${ih}
  done
done
