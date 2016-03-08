#!/bin/sh

SIESTA="$1"
# Find absolute path -------
#
dir=$(dirname $SIESTA)
pushd ${dir}
ABS_EXEC_DIR=$(pwd)
popd
SIESTA=${ABS_EXEC_DIR}/siesta

echo "Running script with SIESTA=$SIESTA"
#
for i in `cat config`; do
  cp -rp ../$i .
  cd $i
  echo $SIESTA < pt2.fdf |tee OUT
  cd ..
done
#

