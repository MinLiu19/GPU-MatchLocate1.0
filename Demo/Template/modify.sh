#!/bin/bash
cat ../catalog.dat|gawk '{print $1}'>temp.list
for event in `cat temp.list`
do
cd $event
sac<<EOF
r *
ch lcalda TRUE
w over
q
EOF
cd ..
done
