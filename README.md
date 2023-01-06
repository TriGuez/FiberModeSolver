# FiberModeSolver
This code simulates the propagation modes of various &amp; complex optical fibers

The solver is based on a finite difference scheme of the following scalar Helmholtz equation : 

$$\big[\Delta_{\perp} + k_{0}^{2}n^{2}(x,y,\omega)\big]\psi(x,y,\omega) = \beta^2(\omega)\psi(x,y,\omega)\ \ \ \ \ \ \ \ \ \  \textbf{(1)}$$
where $\psi$ is the complex electric field, $\beta$ its corresponding propagation constant, $k_{0} = {{2\pi}\over{\lambda}}$ its wavevector and $n$ the refractive index profile of the fiber.

In finite difference domain, assuming that spacing $\Delta x = \Delta y = h$, $\Delta_{\perp}\psi$ can be express as : 
$$\big[\partial_{x}^{2} + \partial_{y}^{2}\big]\psi(x,y) = {1 \over h^2}\big[{f_{p-1,q}+f_{p+1,q}-4f_{p,q}+f_{p,q-1}+f_{p,q+1}}\big]\ \ \ \ \ \ \ \ \ \  \textbf{(2)}$$

Injecting the finite difference scheme into (1) leads to an eigenvalue matrix problem, which can be easily integrate using Matlab function **eigs()**

# Performance comparison

For step index fibers, one can easily compare the solver results with the analytical solution.

We compared effective indexes obtained with this solver with the ones obtained using the semi-analytical improved effective index method [1] for the air silica microstructured fiber described in [2].
As it can be seen in fig.1, results are very similar, leading to a mean squared error below $10^{-8}$. Although the improved effective index method is way faster to estimate the effective index of the fundamental mode of a given air-silica microstructured fiber, this solver 1) gives information about the electric field, and 2) is working for various refractive index profiles, including the ones of hollow core & photonic bandgap fibers.
| ![image](https://user-images.githubusercontent.com/121152666/208872452-8e5f1dd7-0611-402c-805c-b7a6173ab6bf.png)  |
| :--: |
| Fig. 1 : $n_{eff}$ comparison |

# Examples
```Matlab
MMFStandard.m
```
| ![105-125](https://user-images.githubusercontent.com/121152666/208956926-1cf9dd6c-3573-42bd-a34e-04c2bc26257e.gif) |
| :--: |
| First 50 modes of a 105-125, 0.22 NA optical fiber |

```Matlab
PCFStandard.m
```
| ![PCF](https://user-images.githubusercontent.com/121152666/208959172-effad89e-656c-4b68-8438-71d4f5946f5d.gif) |
| :--: |
| First 3 modes of an air-silica microstructured fiber |
```Matlab
ARFStandard.m
```
| ![image](https://user-images.githubusercontent.com/121152666/208961055-6e369851-90c2-4407-8989-519724f215ff.png) |
| :--: |
| Fundamental mode of a generic antiresonnant fiber |

# Functions

##  Refractive Indexes
### - ARFIndex.m
Usage : 
```Matlab
n = ARFIndex(X, Y, FiberParams)
```
Description : 

This function creates the refractive index profile of a fused silica made anti-resonnant fiber, including material optical losses (i.e complex refractive index)

Inputs : 
* X, Y : spatial mesh
* FiberParams : Matlab structure with fields : Dclad : Cladding diameter, NCap : number of capillaries, DextCap : Outter diameter of the capillaries, ECap : Capillaries thickness, lambda : Working wavelength (m)

Output : 
* n : Refractive index map 

### - ASBandgapIndex.m
Usage : 
```Matlab
n = ASBandgapIndex(X, Y, FiberParams)
```

Description : 

This function creates the refractive index profile of an all-solid photonic bandgap fiber (hexagonal lattice of high-index parabolic inclusions), including material optical losses (i.e complex refractive index)

Inputs : 
* X, Y : spatial mesh
* FiberParams : Matlab structure with fields : Pitch : lattice pitch (m), dop : $d\over \Lambda$ parameter, dn : refractive index difference of the high index inclusions, Dclad : Cladding diameter, lambda : working wavelength (m)

Outputs : 
* n : Refractive index map

### - BendedIndex.m
Usage : 
```Matlab
n = BendedIndex(X, Y, RIndexMap, BendingParams)
```

Description : 

This function computes the impact of bending on a given refractive index profile according to [3].

Inputs : 
* X, Y : spatial mesh
* RIndexMap : Input refractive index profile
* BendingParams : Matlab structure with fields : R : Bending radius (m), Angle : Angle between bending direction & the horizontal direction (rad).

Outputs : 
* n : Bended refractive index profile

### - ParabolicIndex.m
Usage :
```Matlab
n = ParabolicIndex(X, Y, FiberParams)
```

Description : 

This function creates the refractive index profile of a parabolic graded index fiber, including material losses (i.e complex refractive index)

Inputs : 
* X, Y : spatial mesh
* FiberParams : Matlab structure with fields : Dcore : Core diameter (m), Dclad : Cladding diameter (m), dn : Refractive index difference of the parabola, lambda : WOrking wavelength (m)

Outputs : 
* n : refractive index profile

### - PCFIndex.m
Usage : 
```Matlab
n = PCFIndex(X, Y, FiberParams)
```
Description : 

This function creates the refractive index profile of an air silica photonic crystal fiber (One hole missing triangular lattice, 5 rows of holes), including material optical losses (i.e complex refractive index)

Inputs : 
* X, Y : spatial mesh
* FiberParams : Matlab structure with fields : Pitch : lattice pitch (m), dop : $d\over \Lambda$ parameter, lambda : Working wavelength (m), Dclad : Cladding diameter

Output : 
* n : Refractive index map 

### - StepIndex.m
Usage : 
```Matlab
n = StepIndex(X, Y, FiberParams)
```
Description : 

This function creates the refractive index profile of step index fiber, including material optical losses (i.e complex refractive index)

Inputs : 
* X, Y : spatial mesh
* FiberParams : Matlab structure with fields : Dcore : Core diameter (m), Dclad : Cladding diameter (m), NA : Core numerical apperture , lambda : Working wavelength (m)

Output : 
* n : Refractive index map 

## Tools

### - ModeArea.m
Usage : 
```Matlab
Aeff = ModeArea(x, y, Field)
```
Description : 

This function computes the effective area of the considered optical field, with : 
$$A_{eff} =  {\Big({\int{\int{\vert E\vert^{2}}} dxdy\Big)^{2}}\over{\int{\int{\vert E\vert^{4}dxdy}}}}\ \ \ \ \ \ \ \ \ \  \textbf{(3)}$$

Inputs : 

* x, y : Transverse coordinates of the simulation (m)
* Field : Matrix of the optical field of the considered mode (must be complex !)

Outputs : 

* Aeff : Effective mode area (m²)

### ModeFieldDiameter.m : 
Usage : 
```Matlab
MFD = ModeFieldArea(x, y, Field, [OPTIONAL]'type', typeValue)
```

Description : 

This function computes the mode field diameter of a given optical mode, using the specified definition (default is 4 $\sigma$)

Inputs : 

* x, y : Transverse coordinates of the simulation (m)
* Field : Matrix of the considered optical field or optical instensity
* [OPTIONAL] type : MFD definition. Options are : '4sigma' (default), '1/e', '1/e2', 'FWHM'

Outputs : 

* MFD : Computed mode field diameter

### ModeSolver.m
Usage :
```Matlab
[neff, LP] = ModeSolver(RIMap, x, y, [OPTIONAL] : 'nModes'(default = 10), nModesValue, 'coreRadius'(default = 5e-6), coreRadiusValue, ...
                     'lambda'(default = 1064e-9, lambdaValue, 'plot'(default = true), plotValue, 'target'(default = 0), targetValue...
                     'IndexContour'(default = true), indexContourValue)
```

Description : 

This script is the main solver, based on finite difference scheme for eq (1)

Inputs : 

* RIMap : Refractive index map of an arbitrary optical fiber
* x, y : Transverse coordinates of the simulation (m)
* [OPTIONAL]
* 'nModes' : Amount of mode to compute. has to be increased if the solver did not converged. Default is 10
* 'coreRadius' : Core radius of the considered fiber, used to evacuate cladding modes. Default is 5e-6
* 'lambda' : Working wavelength. Default is 1064e-9
* 'plot' : Enables the display of the computed modes. Set to false to increase speed, especially for large amount of modes. Default is true
* 'target' : Highest value for the effective index to compute. If set to 0, the eigenvalue solver will look for the highest real part eigenvalue. Default is 0
* 'indexContour' : Enables the display of the optical fiber section over the mode plot. Default is true

Outputs : 

* neff : Effective index of the computed modes
* LP : Corresponding core modes

### SilicaIndex.m
Usage : 
```Matlab
index = SilicaIndex(lambda)
```

Description : 

This function computes the refractive index of fused silica corresponding to the given wavelength, using Sellmeier's formula.

Input : 
* lambda : Working wavelength (m)

Output : 
* Corresponding refractive index

### SilicaLosses.m
Usage : 
```Matlab
alpha = SilicaLosses(lambda)
```

Desctiption : 

This script computes the optical losses in fused silica in regard of the given wavelength, according to the model given in [3].

Input : 

* lambda : Working wavelength (m)

Output : 

* alpha : Optical losses (m $^{-1}$)

# Further improvements
* Perfecly matched layer to evaluate mode leakage losses (WORK IN PROGRESS)
* Fiber bending &#9745;

# References
* [1] K. N. Park et K. S. Lee, *Improved effective-index method for analysis of photonic crystal fibers*, Opt. Lett., vol. 30, nᵒ 9, p. 958, mai 2005, doi: 10.1364/OL.30.000958.

* [2] T. Gottschall et al., *Ultra‐compact tunable fiber laser for coherent anti‐Stokes Raman imaging*, J Raman Spectrosc, vol. 52, nᵒ 9, p. 1561‑1568, sept. 2021, doi: 10.1002/jrs.6171.

* [3] R. T. Schermer and J. H. Cole, "Improved Bend Loss Formula Verified for Optical Fiber by Simulation and Experiment," in IEEE Journal of Quantum Electronics, vol. 43, no. 10, pp. 899-909, Oct. 2007, doi: 10.1109/JQE.2007.903364.

* [4] : Sørensen, S. T. *Deep-blue supercontinuum light sources based on tapered photonic crystal fibers*, PhD thesis, DTU Fotonik ,2013.

# Citation
If this code was helpful, please consider citing it : 


[![DOI](https://zenodo.org/badge/580677271.svg)](https://zenodo.org/badge/latestdoi/580677271)



