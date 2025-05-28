#! /bin/bash

no_particles=2
no_steps=10
Lx=10
Ly=10
h=0.5
v0=10
dt=0.01
c_s=0.01
gamma_index=7

./sph $no_particles $no_steps $Lx $Ly $h $v0 $dt $c_s $gamma_index
