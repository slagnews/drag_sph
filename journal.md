# Weekly progress journal

## Instructions

In this journal you will document your progress of the project, making use of the weekly milestones.

Every week you should 

1. write down **on the day of the lecture** a short plan (bullet list is sufficient) of how you want to 
   reach the weekly milestones. Think about how to distribute work in the group, 
   what pieces of code functionality need to be implemented.
2. write about your progress **until Tuesday, 11:00** before the next lecture with respect to the milestones.
   Substantiate your progress with links to code, pictures or test results. Reflect on the
   relation to your original plan.

We will give feedback on your progress on Tuesday before the following lecture. Consult the 
[grading scheme](https://computationalphysics.quantumtinkerer.tudelft.nl/proj3-grading/) 
for details how the journal enters your grade.

Note that the file format of the journal is *markdown*. This is a flexible and easy method of 
converting text to HTML. 
Documentation of the syntax of markdown can be found 
[here](https://docs.gitlab.com/ee/user/markdown.html#gfm-extends-standard-markdown). 
You will find how to include [links](https://docs.gitlab.com/ee/user/markdown.html#links) and 
[images](https://docs.gitlab.com/ee/user/markdown.html#images) particularly.

## Week 1 - planning the project
(due 21 May 2025, 23:59)
Use this plannning issue to provide a short outline and plan for your final project. Please also include links to existing literature.

We will be implementing a Smoothed Particle Hydrodynamics simulation. We will start of by implementing the basics, with a simple Gaussian or cubic spline interpolation kernel, and checking if mass, momentum and energy are conserved [1].

If all of that works, we can move on to implementing more advanced stuff, such as:
- Calculating the specific heat for an ideal gas equation of state, which can be validated with existing literature values.
- Different equations of state, such as the Cole Equation Of State to model incompressible fluids, or even solids [2]
- Increase the Reynolds number and analyze when turbulence occurs[3] [4], this also allows for easy validation with existing literature

An optional extension could be to implement the computationally heavy parts of the code in C++ with a python wrapper.



[1] Rosswog, S. (2009). Astrophysical smooth particle hydrodynamics. New Astronomy Reviews, 53(4–6), 78–104. https://doi.org/10.1016/j.newar.2009.08.007

[2] Monaghan, J. J. (1992). Smoothed particle hydrodynamics. Annual Review of Astronomy and Astrophysics, 30(1), 543–574. https://doi.org/10.1146/annurev.aa.30.090192.002551

[3] Antuono, M., Marrone, S., Di Mascio, A., & Colagrossi, A. (2021). Smoothed particle hydrodynamics method from a large eddy simulation perspective. Generalization to a quasi-Lagrangian model. Physics of Fluids, 33(1). https://doi.org/10.1063/5.0034568

[4] Jiang, F., Oliveira, M. S., & Sousa, A. C. (2006). SPH simulation of transition to turbulence for planar shear flow subjected to a streamwise magnetic field. Journal of Computational Physics, 217(2), 485–501. https://doi.org/10.1016/j.jcp.2006.01.009

## Week 2
(due 27 May 2025, 11:00)
This week we implemented the basics of the SPH simulation. Starting with a simple non-interacting N body simulation with random motions and periodic boundary conditions. We could re-use quite some code from the first project. Calculating the interpolated densities we start with a simple gaussian kernel, which in 2D is
$$
W(r,h)=\frac{1}{2\pi h^2}\exp{\left(-\frac{1}{2}\frac{r^2}{h^2}\right)}.
$$
Using this, we calculate the density at a grid simply by taking a grid position, calculating its distance to each particle $i$: $r=|\mathbf{r}-\mathbf{r}_{i}|$, and summing over all particles
$$
\rho(\mathrm{r})=\sum_{i}{m_iW(|\mathbf{r}-\mathbf{r}_i|,h)}.
$$
Since we are using periodic boundary conditions, we calculate the distance using
$$
\Delta x=(x-x_i+L/2)\%L-L/2
$$
and
$$
\Delta y=(y-y_i+L/2)\%L-L/2
$$
where $L$ is the size of the simulation space. Now the distance is
$$
r=\sqrt{\Delta x^2+\Delta y^2}.
$$
To check if all went correctly, we calculated the mass first by simply summing over all particle masses
$$
M_1=\sum_i m_i
$$
and secondly by integrating over the density
$$
M_2=\int_0^L\int_0^L \rho(\mathbf{r})\mathrm{d}x\mathrm{d}y.
$$
For 100 particles, with $L=3$, $h=0.5$ and $m_i=1$ for all particles, this of course results in $M_1=100$. The integration over the density resulted in $M_2=99.45$, which is slightly lower. This is to be expected, since the density only takes into account the density contribution of the nearest image of each particle, but since our kernel is Gaussian and extends infinitely, the images that are further away would still have a tiny contribution. It is also partly due to numerical integration over finitely sized volume elements. So this seems to work properly.

### Drag on an object
As a first test of the simulation, we initialized 1000 particles on a 10x10 simulation space with randomly initialized velocities (normally distributed around 0 with an std of `v_random=1`). By giving them an additional `v0=10` velocity in the positive x-direction, the fluid starts to flow to the right.

In the simulation we can add an optional `central_object`: a circular object that does not move and pushes particles to flow around it. This object gives an additional force to each object
$$
F_\mathrm{obj}=\frac{F_\mathrm{max}}{1-\exp{\frac{r-R}{w}}}
$$
with $r$ the distance of the particle to the center of the object, $R$ the object's radius and $w$ the object's boundary width, and $F_\mathrm{max}$ the maximum force of the object. This force is directed radially outwards from the object, such that it repulses particles that come within it's radius.

By doing this, the fluid will slow down due to drag from the central object. The idea of an SPH simulation is that the particles start to behave as a continuous fluid, therefore the drag should be Stokes drag:
$$
F = \frac{1}{2}\rho v^2 Ac_\mathrm{D}
$$
with $\rho$ the density of the fluid around it, $A$ the object's frontal area, $c_\mathrm{D}$ the object's drag coefficient and $v$ the velocity of the object through the fluid. Then the change in velocity is
$$
\frac{\mathrm{d}v}{\mathrm{d}t}=-kv^2, \quad k = \frac{\rho A c_\mathrm{D}}{2m}.
$$
We can solve the differential equation through
$$
\frac{\mathrm{d}{v}}{v^2}=-k\mathrm{d}t
$$
which results in
$$
\frac{1}{v_0}-\frac{1}{v}=-k(t-t_0)
$$
so the velocity as a function of time is
$$
v=\frac{v_0}{1+kv_0t}
$$
where we took $v=v_0$ at $t=0$. Since we are not slowing down the object, but the flow itself, $m$ is the mass of the flow. Furthermore, we can write $\rho=\frac{m}{L^2d}$ with $d$ some imaginary depth (third dimension) of the simulation. Then the frontal area of the object $A$ becomes $A=2Rd$ with $R$ it's radius. So we get
$$
k=\frac{Rc_\mathrm{D}}{L^2}
$$

Below, we show the evolution of the mean of the velocity in the x direction
![](figures/drag.png)
Here we plotted the simulation result in blue, the analytical expectations from Epstein drag in orange (solid line) and a fit from Stokes drag ($F\propto v^2$) in orange (dotted line). As we can see, the results do not quite match. First of all: the analytical result does not slow down nearly enough. This might well be due to the fact that we are using periodic boundary conditions, and therefore particles move past the object again and again, so they effectively encounter it multiple times, which adds up. Looking at the fit we also see that the shape of the curve does not quite seem to match the shape of the data.

What could be happening is that the simulation is not tuned correctly to make the fluid behave as a continuous flow, and instead interacts with the object as a simple N-body code. This would mean that the drag experienced by the object is not Stokes drag, but Epstein drag. Epstein drag is the force on a small object that is about the same size as the mean free path of particles in a gas/fluid (so for life-sized objects Stokes drag is applicable, but for tiny dust particles collisions with individual gas/fluid molecules become important, which is described by Epstein drag). Epstein drag force is given by
$$
F=\frac{4}{3}\rho A v_\mathrm{th}v
$$
such that we get
$$
\frac{\mathrm{d}v}{\mathrm{d}t}=-kv, \quad k \equiv \frac{4\rho A v_\mathrm{th}}{3m}
$$
with $m$ the object's mass, $\rho$ the fluid density, $A$ the object frontal surface area, and $v_\mathrm{th}$ the thermal velocity of particles. We can solve the differential equation by writing
$$
\frac{\mathrm{d}v}{v}=-k\mathrm{d}t
$$
so we get
$$
v=v_0 e^{-kt}
$$
In our simulation, the density $\rho$ is equivalent to a surface density $\sigma$ divided over an imaginary depth $d$, while the surface area is $A=2Rd$, so we get $\rho A=\frac{\sigma}{D}dD=\sigma D$. Now since we are not slowing down the object, but the flow itself, the mass of the object is the total mass of the fluid, which we now can call $m$. Therefore the surface density is $\sigma=m/L^2$ with $L$ the size of the simulation. Using all of this we get
$$
k=\frac{8Rv_\mathrm{th}}{3L^2}
$$
In the image above, we plotted the analytical result as well in green (solid line), which is again not slowing down enough, again due to periodic boundary conditions. But the fit does seem to be quite good, as the overall shape looks similar to the simulation result in blue. Also: the MSE of the Stokes fit is $1.9\times 10^{-2}$ whereas the MSE of the Epstein fit is a lot smaller: $5.8\times 10^{-3}$.

This means that the simulation is creating Epstein drag, meaning that is does not yet corretly simulate a continuous fluid. Possible explanations are that we do not use enough particles, or that we do not include viscosity just yet (as this is important in the Stokes regime). Another thing might be that we are not using the correct smoothing kernel (Gaussian instead of the more accepted splines) and also do not use a large enough smoothing length. These are things that require a bit more research and that we will do next week.


Todo next week:
- [ ] Figure out how to model continuous flow around object (smoothing function, viscosity, number of particles)
- [ ] Implement other equation of state (some simulations of fluids use the ideal gas law with a spline kernel)
- [ ] Look at energy conservation & implement Verlet integration
- [ ] Optimize code (using a finite kernel will allow for great optimization)

Optional:
- [ ] Implement walls & gravity to simulate fluid in a box, an often used example of SPH


## Week 3
(due 3 June 2025, 11:00)


## Reminder final deadline

The deadline for project 3 is **9 June 23:59**. By then, you must have uploaded the presentation slides to the repository, and the repository must contain the latest version of the code.
