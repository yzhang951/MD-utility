#!/bin/bash
date
rm Ni*cfg
mpiexec.mpich -np 32 lmp2020 -in input -partition 32x1> LOG &
