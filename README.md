# FiberModeSolver
This code simulates the propagation modes of various &amp; complex optical fibers

It's based on a finite difference scheme of the following scalar Helmholtz equation : 

$$\big[\Delta_{\perp} + k_{0}^{2}n^{2}(x,y,\omega)\big]\psi(x,y,\omega) = \beta^2(\omega)\psi(x,y,\omega)\ \ \ \ \ \ \ \ \ \  \textbf{(1)}$$
where $\psi$ is the complex electric field, $\beta$ its corresponding propagation constant, $k_{0} = {{2\pi}\over{\lambda}}$ its wavevector et $n$ the refractive index profile of the fiber.

In finite difference domain, assuming that spacing $\Delta x = \Delta y = h$, $\Delta_{\perp}\psi$ can be express as : 
$$\big[\partial_{x}^{2} + \partial_{y}^{2}\big]\psi(x,y) = {1 \over h^2}\big[{f_{p-1,q}+f_{p+1,q}-4f_{p,q}+f_{p,q-1}+f_{p,q+1}}\big]$$

Injecting the finite difference scheme into (1) leads to an eigenvalue matrix problem, which can be easily integrate using Matlab function **eigs()**

## Performance comparison

For step index fibers, one can easily compare the solver results with the analytical solution.

We compared effective indexes obtained with this solver with the ones obtained using the semi-analytical improved effective index method [1] for the air silica microstructured fiber described in [2].
As it can be seen in fig.1, results are very similar, leading to a mean squared error below $10^{-8}$. Although the improved effective index method is way faster to estimate the effective index of the fundamental mode of a given air-silica microstructured fiber, this solver 1) gives information about the electric field, and 2) is working for various refractive index profiles, including the ones of hollow core & photonic bandgap fibers.
| ![image](https://user-images.githubusercontent.com/121152666/208872452-8e5f1dd7-0611-402c-805c-b7a6173ab6bf.png)  |
| :--: |
| Fig. 1 : $n_{eff}$ comparison |

## References
* [1] K. N. Park et K. S. Lee, « Improved effective-index method for analysis of photonic crystal fibers », Opt. Lett., vol. 30, nᵒ 9, p. 958, mai 2005, doi: 10.1364/OL.30.000958.

* [2] T. Gottschall et al., « Ultra‐compact tunable fiber laser for coherent anti‐Stokes Raman imaging », J Raman Spectrosc, vol. 52, nᵒ 9, p. 1561‑1568, sept. 2021, doi: 10.1002/jrs.6171.

