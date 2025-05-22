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

## Week 3
(due 3 June 2025, 11:00)


## Reminder final deadline

The deadline for project 3 is **9 June 23:59**. By then, you must have uploaded the presentation slides to the repository, and the repository must contain the latest version of the code.
