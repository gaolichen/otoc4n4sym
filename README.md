# A command line tool to compute OTOC for N=4 supersymmetric Yang–Mills theory
This is a command line tool to numerically compute Out-of-time-ordered correlator (OTOC) for the [N=4 SYM](https://en.wikipedia.org/wiki/N_%3D_4_supersymmetric_Yang%E2%80%93Mills_theory)
and beta deformed N=4 SYM, for both inifnite N and finite N cases.

All computations are limited to SU(2) sector where trace states are built out of from X and Z. The number of fields in a trace state is called trace length, denoted by integer 
number L. The X field are also called magnon, its number is denoted by an integer M <= L.

For more details of the background, read the paper [One-Loop Non-Planar Anomalous Dimensions in Super Yang-Mills Theory
](https://arxiv.org/abs/2005.14254) and [Regularization dependence of the OTOC. Which Lyapunov spectrum is the physical one?
](https://arxiv.org/abs/1903.09595)

## Features
- Generate trace states of SU(2) sector and Hamiltonian matrix in trace states.
- Generate norm matrix, whose entry at `(i,j)` represents the overlap between trace states `i` and `j`.
- Compute energy spectrum for Hamiltonian matrix
- Compute OTOC


## Build & Run
- Install the [BOOST](https://www.boost.org/) C++ library
- Install [CMake](https://cmake.org/)
- Change `BOOST_ROOT` in CMakeLists.txt in the project root:
```
SET(BOOST_ROOT ${MYDEV}/boost_1_70_0) --> SET(BOOST_ROOT path/toboost/library)
```
- Create a `build` directory under project root
- Run command line `cmake ../` from the `build` directory
- Build the project:
  - In Windows/Visual Studio, double click the .sln file in the `build` directory, then build the solution 
  - In Linux, run command line `make` in the `build` directory
- Run `otoc.exe` from command line to see help & examples

## Examples
Examples are displayed with helps when `otoc.exe` is executed without option.

1. Generate Hamiltonian matrix and trace states for L=4, M=2
```shell
otoc ham -L 4 -M 2
```

2. Generate all Hamiltonian matrix and trace states for L=4
```
otoc ham -L 4
```

3. Generate norm matrix for L=4, M=2
```
otoc norm -L 4 -M 2
```

4. Generate all norm matrix for L=4
```
otoc norm -L 4
```

5. Compute energy spectrum for L=4, M=2, N=17, beta=0.9
```
otoc spec -L 4 -M 2 -N 17 -b 0.9
```

6. Generate matrix for operator counting traces containing X for L=4
```
otoc operator -o NTrX -L 4
```

7. Compute OTOC for W=Tr(Xd/dZ), V=Tr(Zd/dX), L=4, N=5, beta=0.9, beta of temperature=0.5, and time run from 0.0 to 5.0 with step size 0.1. Note that
to run the example one needs to run `otoc norm - L 4` first to generate all norm matrices for L = 4 and then copy the generated files to otocdata folder
```
otoc otoc -W XZ -V ZX -L 4 -N 5 -b 0.9 -bt 0.5 -tmin 0.0 -tmax 5.0 -tstep 50
```

8. Compute OTOC for W=Tr(XZ d/dX d/dZ), V=Tr(ZX d/dX d/dZ), L=4, M=2, beta=0.9, beta of temperature=0.5, and N=Infinity. As -n option is specified,
it will compute norm matrix rather than load it from local files.
```
otoc otoc -W XZXZ -V XZXZ -L 4 -M 2 -b 0.9 -bt 0.5 -tmin 0.0 -tmax 5.0 -tstep 50 -n
```

9. Compute C(t) for W=Tr(XX d/dX d/dZ), V=Tr(ZZ d/dZ d/dX), L=4, N=5. Zero energy states are excluded.
```
otoc otoc -W XXXZ -V ZZZX -L 4 -N 5 -b 0.9 -bt 0.5 -tmin 0.0 -tmax 5.0 -tstep 50 -nozero
```

10. Compute C(t) for W=Tr(XX d/dX d/dZ), V=Tr(ZZ d/dZ d/dX), L=4, N=5. With regularization parameters alpha=0.5 and sigma=0.25.
```
otoc otoc -W XXXZ -V ZZZX -alpha 0.5 -sigma 0.25 -L 4 -N 5 -b 0.9 -bt 0.5 -tmin 0.0 -tmax 5.0 -tstep 50
```

11. Compute OTOC for composite operators. Currently, addition (+), substraction(-), multiplication(*), and exponentiate (exp), and normal ordered (~) are
supported. No division (/) and space is allowd in the expression. Complex number should write in the form (real,imag).
```
otoc otoc -W -0.2*exp(0.2*XXXZ+0.3*XZ)+(0.3,-0.5) -V ~(-XZXZ+0.5*ZZZX-0.2*ZX) -L 4 -N 5 -b 0.9 -bt 0.5 -tmin 0.0 -tmax 5.0 -tstep 50
```
