# FiberModeSolver
This code simulates the propagation modes of various &amp; complex optical fibers

It's based on a finite difference scheme of the following Helmholtz equation : 

$$\big[\Delta_{\perp} + k_{0}^{2}n^{2}(x,y,\omega)\big]\psi(x,y,\omega) = \beta^2(\omega)\psi(x,y,\omega)\ \ \ \ \ \ \ \ \ \  \textbf{(1)}$$
where $\psi$ is the complex electric field, $\beta$ its corresponding propagation constant, $k_{0} = {{2\pi}\over{\lambda}}$ its wavevector et $n$ the refractive index profile of the fiber.

In finite difference domain, assuming that spacing $\Delta x = \Delta y = h$, $\Delta_{\perp}\psi$ can be express as : 
$$\big[\partial_{x}^{2} + \partial_{y}^{2}\big]\psi(x,y) = {1 \over h^2}\big[{f_{p-1,q}+f_{p+1,q}-4f_{p,q}+f_{p,q-1}+f_{p,q+1}}\big]$$

Injecting the finite difference scheme into (1) leads to a 

![image](https://user-images.githubusercontent.com/121152666/208872452-8e5f1dd7-0611-402c-805c-b7a6173ab6bf.png)
