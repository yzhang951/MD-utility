#!/bin/bash
date
mpiexec -np 60 lmp -in input > LOG &
