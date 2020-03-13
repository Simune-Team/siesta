#!/bin/bash

# Stop at first error
set -e

# Init
siesta_bindir="/siesta/bin"
siesta_cfgdir="/siesta"
siesta_srcdir="/siesta/trunk"
siesta_version="4.0.2"
vendor="gnu"

# Update the source tree
cd "${siesta_srcdir}"
git checkout master
git clean -fd
git clean -fdX
git pull --prune
git checkout "v${siesta_version}"

# Install selected flavors of SIESTA
mkdir -p "${siesta_bindir}"
for flavor in serial mpi; do
  siesta_builddir="${siesta_srcdir}/Obj/${flavor}"
  mkdir -p "${siesta_builddir}"
  cd "${siesta_builddir}"
  /bin/bash "${siesta_srcdir}/Src/obj_setup.sh"
  cp "${siesta_cfgdir}/arch.make-siesta-${siesta_version}-${vendor}-${flavor}" \
    "arch.make"
  make clean
  make
  cp "siesta" "${siesta_bindir}/siesta-${siesta_version}-${flavor}"
done
