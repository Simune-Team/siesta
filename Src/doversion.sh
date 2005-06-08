#/bin/sh
#
# Simple filter to insert the arch version string
# For some strange reason I could not manage
# to do it directly in the Makefile
#
# This script acts like a filter, replacing the
# string SIESTA_VERSION by the contents of the
# file "version.info", which is automatically generated
# by an arch pre-commit hook.
#
read version < ../version.info
sed  "s'SIESTA_VERSION'${version}'g" 
