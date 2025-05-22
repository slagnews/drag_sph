# Project 3
The central idea of Smoothed Particle Hydrodynamics (SPH) simulations is to approximate continuous materials by a set of particles, of which the properties are interpolated over the space in between them. This is done using a kernel function $W(\mathbf{r}-\mathbf{r}_\mathrm{pt},h)$ where $\mathbf{r}$ is the position at which the property is being evaluated, and $\mathbf{r}_\mathrm{pt}$ the position of the particle of which the contribution is being assessed, and $h$ is a scale factor which is a measure for how spread out the properties are.

In an abstract sense, any property $A(\mathbf{r})$ (this could be the density), can be approximated by sampling at different positions $mathbf{r'}$ and spreading these out using such a kernel function, such that we end up with an interpolated approximation
$$
A_\mathrm{int}(\mathbf{r})=\int_\mathcal{V}A(\mathbf{r}')W(\mathbf{r}-\mathbf{r}',h)\mathrm{d}\mathbf{r}'
$$

The acceleration of a particle $i$ is given by
$$
\mathbf{a}_i=-\sum_{j}m_j \left(\frac{P(\mathbf{r}_i)}{\rho(\mathbf{r}_i)^2}+\frac{P(\mathbf{r}_j)}{\rho(\mathbf{r}_j)^2}\right)\nabla W_{ij}
$$
(see Monaghan eq. 3.8). with $\nabla W_{ij}$ the gradient of the kernel function of $j$ at the position of $i$, so
$$
\nabla W_{ij}=\frac{\mathbf{r}_i-\mathbf{r}_j}{|\mathbf{r}_i-\mathbf{r}_j|}\frac{\mathrm{d}}{\mathrm{d}r}W(r,h)|_{r=|\mathbf{r}_i-\mathbf{r}_j|}
$$