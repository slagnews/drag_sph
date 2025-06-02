#! /bin/bash

res=0.15
no_steps=1000
Lx=5
Ly=5
rho0=1000
kappa=1.2
v0=10
dt=0.01
c_s=0.01
gamma_index=7

./sph $res $no_steps $Lx $Ly $rho0 $kappa $v0 $dt $c_s $gamma_index
