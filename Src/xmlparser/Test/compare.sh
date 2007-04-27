#!/bin/bash

reference=h2o.xml;

number_of_tests=16;
i=1;

while [ $i -le $number_of_tests ]
do
	echo "Testing file: " $i
	../test-xml $reference h2o-$i.xml
        echo "-------------------------";
	let i=i+1;
done
