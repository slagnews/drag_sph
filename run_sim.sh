#! /bin/bash

res=0.1
no_steps=1000
Lx=5.0
Ly=5.0
rho0=1.0
kappa=1.2
v0=1.0
dt=0.000005
c_s=200
beta=0.07
central_radius=0.5
boundary_width=0.02
max_force=0
inflow_factor=0.1
outflow_factor=0.1
kinematic_viscosity=1
max_beta=1.5

./src/sph $res $no_steps $Lx $Ly $rho0 $kappa $v0 $dt $c_s $beta $central_radius $boundary_width $max_force $inflow_factor $outflow_factor $kinematic_viscosity $max_beta