# hierarchicaldd
A series of experiments for the hierarchical domain decomposition algorithm on various platforms.

## Abstract

This description contains the information required to run computational experiments described by the
SC17 paper titled "A Hierarchical Domain Decomposition Preconditioner for the Poisson Equation Optimized for Heterogeneous Platforms". It contains a set of instructions for acquiring, installing, and running the experiments, as a part of ACM reproducibility initiative.


## Description

### Check-list (artifact meta information)

- **Algorithm:** Hierarchical domain decomposition for solving Poisson problems
- **Program:** CUDA and x86/KNL SIMD optimized C++ code
- **Compilation:** x86, KNL C++ code: Intel Compiler (required, v17.0.0), CUDA code: NVIDA nvcc (tested, v5 and v8) with Intel icc (v 17.0.0)
- **Binary:** CUDA, x86, and KNL executables
- **Data Set:** Procedurally generated
- **Runtime Environment:** Ubuntu 14.04, with CUDA toolkit and Intel Parallel Studio XE 2017 (v8 and v2017, respectively)
- **Hardware:** For CUDA experiments, any CUDA enabled GPU with compute 3.5 capability (tested Titan X series). For KNL experiments, a Knights Landing platform with High Bandwidth Memory support. For CPU experiments, any Haswell or newer Intel CPU.
- **Output:** Convergence data represented as geometric visualization and numerical timing.
- **Experiment Work-flow:** Download software and dependencies, build, generate data-sets, run solver experiments, view results with included visualization tool.
- **Publicly Available:** Yes.

### How software can be obtained (if available}

The CPU, GPU, and KNL versions of the experiments can be found in the hierarchicaldd project hosted on Github. Data-sets are generated procedurally from analytic domains. The software comprises code and CMake build files, and is provided under a BSD-3 license.

### Hardware dependencies

For full reproducibility, a system with multiple GPUs, each capable of at least compute v3.5, a Intel Knights Landing platform, and a CPU platform with at least Intel's Haswell architecture or later should be available. 

### Software dependencies

- CMake (>3.0)
- Intel Compiler (2017)
- Eigen (>3.3.1) - http://eigen.tuxfamily.org/
- GLUT
- OpenGL
- Memkind (>1.5.0) - http://memkind.github.io/memkind/
  * (Optional - For Knights Landing Platform)
- CUDA - (>5.0)
  * (Optional - For GPU Implementation)


## Installation

### Download

Retrieve core experimental code from the hierarchicaldd project:

```
git clone https://github.com/uwgraphics/hierarchicaldd
```

Create a build path and execute CMake for the project:

```
cd hierarchicaldd
mkdir build
cd build
cmake ..
```

You need to manually set paths for Intel's icc compiler, the Eigen library, and the memkind library by running ccmake after the initial configuration. Additionally, the GPU and KNL experiments can be enabled through options at this point.

```
ccmake .
```

Finally, the experiments can be built by executing **make** from the build directory.

```
make
```

Once built, the binaries can be found at the following paths:

- tools/opengl\_3d - Visualization tool for 2D convergence result.
- demos/SPGrid\_Create\_Simple_Domain - A project for creating an simple analytic domain with a sphere in the middle and initialize the right hand side with a Kronecker delta function.
- demos/Domain\_Decomposition\_CPU_Only - The CPU implementation of the domain decomposition algorithm proposed by Lui et al\cite{}.
- demos/Domain\_Decomposition\_CPU\_Only\_Two\_Levels - The CPU implementation of the two level domain decomposition algorithm.
- demos/Domain\_Decomposition\_KNL - The KNL implementation of the domain decomposition algorithm.
- demos/Domain\_Decomposition\_CUDA\_Matrix\_Free - The GPU(CUDA) implementation of the domain decomposition algorithm.

## Experiment work-flow

### SPGrid\_Create\_Simple\_Domain:

#### Description
Inside SPGrid_Create\_Simple\_Domain, there are two projects: Create\_3D\_Domain and Create\_2D\_Domain, each cooresponding to create domain of specified dimension.

#### Usage
  - 3D: `./create_simple_domain_nocona [-size <n n n>] [-position_delta <n n n>] [-d <Directory>]`
  - 2D: `./create_simple_domain_nocona [-size <n n>] [-position_delta <n n>] [-d <Directory>]`
  
#### Example
  - 3D: `./create_simple_domain_nocona -size 256 256 256 -position_delta 64 64 64 -d output_256`
  - 2D: `./create_simple_domain_nocona -size 256 256 -position_delta 64 64 -d output_256`
   
### Domain\_Decomposition\_CPU\_Only:

#### Description
Inside Domain\_Decomposition\_CPU\_Only, there are two executable projects: Test\_Bench\_3D and Test\_Bench\_2D, each cooresponding to solve domain of specified dimension.

#### Usage:
- 3D: `./dd_pcg_nocona [-size <n n n>] [-sub_size <n n n>] [-threads <n>] [-vlevels <n>] [-sigma_levels <n>] [-sigma_smoothing <n>] [-sigma_mg <n>] [-interior_smoothing <n>] [-boundary_smoothing <n>] [-nVcycles <n>] [-d <Directory>] -mode no_convertion`
- 2D: `./dd_pcg_nocona [-size <n n>] [-sub_size <n n>] [-threads <n>] [-vlevels <n>] [-sigma_levels <n>] [-sigma_smoothing <n>] [-sigma_mg <n>] [-interior_smoothing <n>] [-boundary_smoothing <n>] [-nVcycles <n>] [-d <Directory>] -mode no_convertion` 

#### Parameters
   * -d : domain file directory
   * -size : size of the domain + 1
   * -sub_size : size of the subdomain
   * -threads : number of threads is being used
   * -vlevels : subdomain solver, v cycle levels
   * -nVcycles : subdomain solver, number of v cycles per subdomain solve
   * -interior_smoothing : subdomain solver, interior smoothing iteration for subdomain v cycele
   * -boundary_smoothing : subdomain solver, boundary smoothing iteration for subdomain v cycele
   * -sigma_levels : interface solver, v cycle levels
   * -sigma_smoothing : interface solver, smoothing iterations for interface v cycle
   * -sigma_mg : interface solver, number of v cycles per interface solve

#### Example
- 3D: `./dd_pcg_nocona -d ../../SPGrid_Create_Simple_Domain/Create_3D_Domain/output_256/ -mode no_convertion -size 257 257 257 -sub_size 128 128 128 -vlevels 4 -sigma_levels 5 -nVcycles 5 -interior_smoothing 1 -boundary_smoothing 5 -sigma_smoothing 1 -sigma_mg 1 -threads 6`
- 2D: `./dd_pcg_nocona -d ../../SPGrid_Create_Simple_Domain/Create_2D_Domain/output_256/ -mode no_convertion -size 257 257 -sub_size 128 128 -vlevels 4 -sigma_levels 5 -nVcycles 5 -interior_smoothing 1 -boundary_smoothing 5 -sigma_smoothing 1 -sigma_mg 1 -threads 6`

### Domain\_Decomposition\_CPU\_Only\_Two\_Levels

#### Description
Inside Domain\_Decomposition\_CPU\_Only\_Two\_Levels, there are two executable projects: Test\_Bench\_3D and Test\_Bench\_2D, each cooresponding to solve domain of specified dimension. This runs very similar to CPU only version, except with an additional flag that indicates the use of the two level domain decomposition algorithm instead of the original algorithm.

#### Usage:
- 3D: `./dd_pcg_nocona [-size <n n n>] [-sub_size <n n n>] [-threads <n>] [-vlevels <n>] [-sigma_levels <n>] [-sigma_smoothing <n>] [-sigma_mg <n>] [-interior_smoothing <n>] [-boundary_smoothing <n>] [-nVcycles <n>] [-d <Directory>] -mode no_convertion` 
- 2D: `./dd_pcg_nocona [-size <n n>] [-sub_size <n n>] [-threads <n>] [-vlevels <n>] [-sigma_levels <n>] [-sigma_smoothing <n>] [-sigma_mg <n>] [-interior_smoothing <n>] [-boundary_smoothing <n>] [-nVcycles <n>] [-d <Directory>] -mode no_convertion [-multilevel <yes or no>]`

#### Parameters:
   * -d : domain file directory
   * -size : size of the domain + 1
   * -sub_size : size of the subdomain
   * -threads : number of threads is being used
   * -vlevels : subdomain solver, v cycle levels
   * -nVcycles : subdomain solver, number of v cycles per subdomain solve
   * -interior_smoothing : subdomain solver, interior smoothing iteration for subdomain v cycele
   * -boundary_smoothing : subdomain solver, boundary smoothing iteration for subdomain v cycele
   * -sigma_levels : interface solver, v cycle levels
   * -sigma_smoothing : interface solver, smoothing iterations for interface v cycle
   * -sigma_mg : interface solver, number of v cycles per interface solve
   * -multilevel : indicates the use of the two level domain decomposition algorithm
   
#### Example
- 3D: `./dd_pcg_nocona -d ../../SPGrid_Create_Simple_Domain/Create_3D_Domain/output_256/ -mode no_convertion -size 257 257 257 -sub_size 128 128 128 -vlevels 4 -sigma_levels 5 -nVcycles 5 -interior_smoothing 1 -boundary_smoothing 5 -sigma_smoothing 1 -sigma_mg 1 -threads 6 -multilevel yes`
- 2D: `./dd_pcg_nocona -d ../../SPGrid_Create_Simple_Domain/Create_2D_Domain/output_256/ -mode no_convertion -size 257 257 -sub_size 128 128 -vlevels 4 -sigma_levels 5 -nVcycles 5 -interior_smoothing 1 -boundary_smoothing 5 -sigma_smoothing 1 -sigma_mg 1 -threads 6 -multilevel yes`

### Domain\_Decomposition\_KNL

#### Description
The KNL implementation takes the same parameters as Domain\_Decomposition\_CPU_Only. But given the specialized optimization required for the Xeon Phi platform, only 3D is supported for this project.

### Domain\_Decomposition\_CUDA\_Matrix\_Free

#### Description
The CUDA implementation takes the same parameters as Domain\_Decomposition\_CPU\_Only. But given the specialized optimization required for the GPU platform, only 3D is supported for this project. nVcycles,interior\_smoothing, and boundary\_smoothing, these three parameters are hard coded in the CUDA implementation.

### opengl\_3d:

#### Description
The visualization tool for viewing the 2D convergence result as height field

#### Example
`opengl_3d ../../demos/Domain_Decomposition_CPU_Only/Test_Bench_2D/output_residual`


## Evaluation and expected result

The expected results include geometric output representing the overall convergence of the solvers, along with numerical timing information displayed on the console while the experiments are running.


## Notes

The intent for the implementation is to optimize the domain decomposition preconditioner call as much as possible. But many CPU implementation are still sub-optimal, other CG kernels, for instance.

For up-to-date instructions, please see the Github project page (https://github.com/uwgraphics/hierarchicaldd).
