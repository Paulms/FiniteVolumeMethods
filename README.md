# Finite Volume Methods

Finite Volume Methods for Hyperbolic Problems. These PDEs are of the form

<a href="https://www.codecogs.com/eqnedit.php?latex=u_{t}&plus;f(u)_{x}&=0,\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{&plus;}\\u(x,0)&=u_{0}(x),\qquad\forall&space;x\in\mathbb{R}^{n}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{t}&plus;f(u)_{x}&=0,\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{&plus;}\\u(x,0)&=u_{0}(x),\qquad\forall&space;x\in\mathbb{R}^{n}" title="u_{t}+f(u)_{x}&=0,\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{+}\\u(x,0)&=u_{0}(x),\qquad\forall x\in\mathbb{R}^{n}" /></a>

We also consider degenerate convection-diffusion systems of the form:

<a href="https://www.codecogs.com/eqnedit.php?latex=u_{t}&plus;f(u)_{x}&=(B(u)u_{x})_{x},\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{&plus;}\\u(x,0)&=u_{0}(x),\qquad\forall&space;x\in\mathbb{R}^{n}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{t}&plus;f(u)_{x}&=(B(u)u_{x})_{x},\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{&plus;}\\u(x,0)&=u_{0}(x),\qquad\forall&space;x\in\mathbb{R}^{n}" title="u_{t}+f(u)_{x}&=(B(u)u_{x})_{x},\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{+}\\u(x,0)&=u_{0}(x),\qquad\forall x\in\mathbb{R}^{n}" /></a>

Solutions follow a conservative finite diference (finite volume) pattern. This method updates point values (cell averages) of the solution **u** and has the general form

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{d}{du}u_{i}(t)=-\frac{1}{\Delta_{i}x}(F_{i&plus;1/2}(t)-F_{i-1/2}(t))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{d}{du}u_{i}(t)=-\frac{1}{\Delta_{i}x}(F_{i&plus;1/2}(t)-F_{i-1/2}(t))" title="\frac{d}{du}u_{i}(t)=-\frac{1}{\Delta_{i}x}(F_{i+1/2}(t)-F_{i-1/2}(t))" /></a>

Where the numerical flux <a href="https://www.codecogs.com/eqnedit.php?latex=F_{i&plus;1/2}(t)&space;=&space;F(u_{i}(t),u_{i&plus;1}(t)))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F_{i&plus;1/2}(t)&space;=&space;F(u_{i}(t),u_{i&plus;1}(t)))" title="F_{i+1/2}(t) = F(u_{i}(t),u_{i+1}(t)))" /></a> is an approximate solution of the Riemann problem at the cell interface ($x_{i+1/2}$). 

An extra term **P** similar to **F** could be added to account for the Diffusion in the second case.

The time integration of the semi discrete form is performed with Runge Kutta methods.

## Features
### Mesh: 
At the momento only Cartesian 1D uniform mesh available, using `FVMesh(N,a,b,boundary)` command. Where

`N` = Number of cells

`xinit,xend` = start and end coordinates.

`bdtype` = boundary type (ZERO_FLUX, PERIODIC)

* Problem types: System of Conservation Laws with diffusion term (`CLS1DDiuffsionProblem`).

### Algorithms

* High-Resolution Central Schemes:

Kurganov, Tadmor, *New High-Resolution Central Schemes for Nonlinear Conservation Laws and Convection–Diffusion Equations*, Journal of Computational Physics, Vol 160, issue 1, 1 May 2000, Pages 241-282

* Entropy Stable Schemes for degenerate convection-diffusion equations

Jerez, C. Pares. *Entropy stable schemes for degenerate convection-difusion equations*. 2017. Society for Industrial and Applied Mathematics. SIAM. Vol. 55. No. 1. pp. 240-264

### Time integration methods:

At the moment available methods are: Forward Euler (`FORWARD_EULER`), Strong Stability Preserving Runge Kutta 2 (`SSPRK22`), `SSPRK33`.

## Example


```fortran
  ! Setup Mesh
  USE FVtypes
  USE EC_scheme
  USE FV_Solve
  INTEGER|                        :: N        ! Number of cells
  REAL(kind = dp)                 :: xinit, xend 
  INTEGER                         :: bdtype
  type(Uniform1DMesh)             :: mesh
  N = 5; xinit = 0.0; xend = 10.0; bdtype = PERIODIC
  call mesh%Initialize(N, xinit, xend, bdtype)

  ! Setup problem
  INTEGER                         :: M    ! number of equations in system
  REAL(kind = dp)                 :: Tend ! Final time
  REAL(kind = dp), ALLOCATABLE    :: uinit(:,:) ! Initial condition
  type(CLS1DDiuffsionProblem)     :: prob
  M = 4; Tend = 0.2;
  ALLOCATE(uinit(N,M))
  uinit = 0.0_dp;
  ! Flux, JcF, BB are functions (see test1.f90 for more information)
  CALL prob%Initialize(mesh, uinit, M, Tend, Flux, JacF, BB)

  ! Choose algorithm and solve problem (some don't require initialization)
  type(ESJP1DAlgorithm)       :: ESPJAlg
  REAL(kind = dp)             :: CFL
  CFL = 0.2
  ! NFlux, KKN are functions (see test1.f90 for more information)
  CALL ESPJAlg%Initialize(Nflux, KKN, 0.0_dp)  
  CALL solve(prob, SSPRK22, ESPJAlg, CFL)
```

# Disclamer
** Modules developed for personal use, some of them have not been tested enough !!!**
