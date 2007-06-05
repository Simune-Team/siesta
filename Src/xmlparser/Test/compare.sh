#!/bin/bash

TESTER=../test-xml
OPTIONS=$@
#
reference=h2o.xml;
echo "---> Using options: $OPTIONS"

number_of_tests=16;
i=1;

while [ $i -le $number_of_tests ]
do
	echo "Testing file: " $i
	$TESTER $OPTIONS $reference h2o-$i.xml
        echo "-------------------------";
	let i=i+1;
done
