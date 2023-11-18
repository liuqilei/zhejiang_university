#!/bin/bash

echo 0 | gmx_mpi trjconv -f pro.trr -s pro.tpr -pbc atom -sep -b 5000 -e 50000 -skip 10 -o system.gro
python3 count.py 983 33 >> CN_condition.dat
rm *.gro
