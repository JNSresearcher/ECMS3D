
---



### Basic mathematical relationships
### 1 Air region
<a id="eq1">**Air region equations:**</a>  
&emsp; Partial differential equations are solved in the time domain for the vector magnetic potential **$\vec A$**

<table><tr> 
<td>
$$\nabla^2 \vec{A} = -\mu_r\mu_0 \vec{J}_s$$ 
</td>
</tr></table> 

where  **$\mu_0=4\pi\cdot 10^{-7}$** is magnetic constant , **$\mu_r=1$** relative magnetic permeability of air (hereinafter simply **$\mu$**),  **$\vec{J}_c$**  is density of external current source.  

<a id="eq1">**Calculation of current density in coils**</a>  
&emsp; The source of the field is the current carrying conductors. The input data define only the current amplitude as a function of time and the geometric configuration of the conductor. The direction of current in conductors and its density is calculated using the following algorithm.  
First, some scalar field is calculated in the coil region using the Laplace equation: **$\nabla^2 V =0$** with Dirichlet boundary conditions at the ends of the current carrying conductor $V_+=+V_0, V_-=-V_0$,  where $+V_0, -V_0$  is the given constants, and the Neumann conditions on the lateral surfaces  $\frac{ \partial V}{\partial n}=0$, &emsp; then the gradient is calculated: $\vec{G} =-\nabla V$, after which the direction cosines are calculated for each cell in the conductor region: $cos_x=\frac{G_x}{|G|}, cos_y=\frac{G_y}{|G|}, cos_z=\frac{G_z}{|G|}$. These values ​​are then multiplied by the time function value obtained from the input data. The algorithm does not automatically calculate the cross-sectional area along the entire length of the conductor. Only the area at the ends of the conductor is calculated automatically. The final current density is obtained by dividing the calculated current by this area.

### 2 Domains of magnetostatic
<a id="eq1">**Equations of magnetostatic region (non-conducting region) :**</a>  
&emsp; Partial differential equations are solved in the time domain for the vector magnetic potential **$\vec A$**  
<table><tr> 
<td>

$$\nabla\times \left( \nabla \times  \vec{A} \right)  = \mu_0 \vec{J}_s$$ 
</td>
</tr></table>  

&emsp; In  initial data on define of the geometric configuration of the regions where the relative magnetic permeability **$\mu_r \neq 1$**.  Note, that for these regions $\vec{J}_s=0$. At the boundaries of regions with different relative magnetic permeabilities **$\mu_1 \neq \mu_2$**, equations for the tangential components of the vector magnetic potential are applied  

<table><tr> 
<td>

$$\frac{1}{\mu_1}\left( \nabla \times  \vec{A_1} \right)  =  \frac{1}{\mu_2}\left( \nabla \times  \vec{A_2} \right) $$ 
</td>
</tr></table>

If the current source borders on a region with relative magnetic permeability **$\mu \neq 1$**, then the equations for the tangential components of the vector magnetic potential are also applied:  

<table><tr> 
<td>

$$\left( \nabla \times  \vec{A_1} \right)  -  \left( \nabla \times  \vec{A_2} \right) =\mu \mu_0 \left( \vec{J}_s \times \vec{n} \right)$$ 

</td>
</tr></table>

where **$\vec{n}$** is the normal to the boundary with the domain, where **$\mu \neq 1$**  
**Note**  This boundary condition is currently under development, so in this version there must be an air gap between external current sources and regions with $\mu \neq 1$.


### 3 Domains of eddy current
<a id="eq2">**Equations of conducting region:**</a>  
&emsp;  Partial differential equations are solved in the time domain for the vector magnetic potential $\vec A$ and the scalar electric potential $U$. 

<table><tr> 
<td>

$$\nabla^2 \vec{A}-\mu_0\sigma \left(\frac{\partial \vec{A}}{\partial t} +\nabla U + (\vec{V}_e \cdot\nabla) \vec{A}  \right)= 0$$ 
  
$$\nabla^2 U + \nabla \cdot \left( \frac{\partial \vec{A}}{\partial t} + (\vec{V}_e\cdot\nabla) \vec{A} \right) = 0$$ 

</td>
</tr></table>

where **$\sigma$**  is conductivity, **$\vec{V}_e$** is velocity of the conducting region in Euler coordinates. The input data can use the velocity of the coils (without rotation) in Lagrangian coordinates **$\vec{V}_s$** (but without crossing the regions where **$\mu \neq 1$** ). If the coils velocity **$\vec{V}_s$**  is specified, then the conducting region is automatically assigned the velocity in Euler coordinates **$\vec{V}_e = -\vec{V}_s$**.  
&emsp; In numerical calculations for a non-conducting region, to determine the **$x,y,z-$** coordinates of the location of external current sources in the coils at time step $i+1$, equations of the form **$x_{i+1}=x_i + Vs_x\cdot dt$** are used, where **$Vs_x$** is the speed of movement of the coil along the coordinate **$x$**, **$dt$** - calculated time step. Similarly, for the remaining coordinates, the speeds of mechanical movement of the coils must be specified: **$Vs_y, Vs_z$**.  
&emsp;  Regions with different conductivities and different velocities can be used, but these regions must not touch (additional boundary conditions are needed here, which are not yet available)

#### 3.1 Equation for calculating eddy current
<table><tr> 
<td>

$$\vec{J_e} = \sigma \left( \frac{\partial \vec{A}}{\partial t} +\nabla U + (\vec{V}_e\cdot\nabla) \vec{A} \right)$$
</td>
</tr></table>
         
&emsp;For clean the eddy currents from divergence (it is required that $\nabla \cdot \vec{J_e}=0$ ). To eliminate divergence, the Hodge-Helmholtz projection method was used:
<table><tr> 
<td>

$$\nabla^2 \phi = \nabla \cdot \vec{J_e}^n;  \vec{J_e}^{n+1} = \vec{J_e}^n - \nabla\phi$$

</td>
</tr></table>

&emsp; To derive these equations, I used the information given in the **Bibliography** [[1]](#b1), [[2]](#b2), [[3]](#b3).   

#### 3.2 Matrix form of equations for calculations:

$$
\begin{aligned}
\begin{bmatrix}
\nabla^2-\mu_0\sigma \frac{\partial}{\partial t} &  &  & -\mu_0\sigma\frac{\partial}{\partial x}  \\
     & \nabla^2 -\mu_0\sigma \frac{\partial}{\partial t}  &  & -\mu_0\sigma\frac{\partial}{\partial y}  \\
  &   & \nabla^2 -\mu_0\sigma \frac{\partial}{\partial t}  & -\mu_0\sigma\frac{\partial}{\partial z}  \\
\frac{\partial}{\partial x} \frac{\partial}{\partial t} & \frac{\partial}{\partial y} \frac{\partial}{\partial t} & \frac{\partial}{\partial z} \frac{\partial}{\partial t} & \nabla^2
\end{bmatrix} \cdot
\begin{bmatrix}
A_x  \\
A_y  \\
A_z \\
U 
\end{bmatrix} =
\begin{bmatrix} 
\mu_0 J_{sx}+\mu_0\sigma\vec{V}_e\cdot\nabla A_x  \\
\mu_0 J_{sy}+\mu_0\sigma\vec{V}_e\cdot\nabla A_y   \\
\mu_0 J_{sz}+\mu_0\sigma\vec{V}_e\cdot\nabla A_z  \\
-\nabla \cdot (\vec{V}_e\cdot\nabla) \vec{A}
\end{bmatrix} 
\end{aligned}
$$

&emsp; The components of the external source **$Js_x, Js_y, Js_z$** are located only in the regions where  **$\sigma=0$**.
&emsp;In this matrix equation, the velocity-dependent components are placed on the right-hand side. The values ​​of the vector potential in the vector of the right sides are determined from the results obtained in the previous calculation step.    
&emsp; Boundary conditions: open boundaries for magnetic vector potential.   For the electric scalar potential, zero Neumann conditions for the normal component at the boundary of the conducting region. Also at the boundary of the conducting region, zero normal components for eddy currents.  
&emsp;Using finite difference approximation, the equations are reduced to a sparse system of algebraic equations. 

### 4 Solvers 
* The direct method is used to solve the governing equations (as well as for test calculations). This method is also used for research purposes to evaluate other, usually iterative, methods.
* To solve auxiliary problems related to calculating the direction of current in coils and eliminating the divergence of eddy currents, a simple iterative method is used.  


***
### Bibliography
<a id="b1">[1]</a>  S. Yamamura, “Theory of linear induction Motors”, John Whiley & Sons, 1972.
ISBN 9780470970904 URL https://www.amazon.com/Theory-Linear-Induction-Motors-Yamamura/dp/0470970901  
<a id="b2">[2]</a>  F.F. Mende “Consideration and the Refinement of Some Laws and Concepts of Classical Electrodynamics and New Ideas in Modern Electrodynamics”,  International Journal of Physics, 2014, Vol. 2, No. 6, 231-263. URL https://pubs.sciepub.com/ijp/2/6/8  
<a id="b3">[3]</a> J. Brackbill, D. Barnes, The effect of nonzero $\nabla \cdot B$ on the numerical solution of the magnetohydrodynamic equations, Journal of Computational Physics 35 (3) (1980) 426–430. doi:10.1016/0021-9991(80)90079-0. URL http://dx.doi.org/10.1016/0021-9991(80)90079-0  

 ***
Autor <a href="mailto:JNSresearcher@gmail.com">J.Sochor</a>

