## Compilation

The dependent libraries are packed into a conda environment, you can create it by

```bash
conda env create -f environment.yml
conda activate msf
```

Then we use cmake to configure and build the project:

```
mkdir build
cd build
cmake .. -DUSE_CONDA=ON -DCMAKE_BUILD_TYPE=Release
make -j8
```

The compilation has been tested on Linux with g++ 9.3.0.



## Usage

#### Optimize asymptotic directional stiffness

```bash
minsurf -j tailads -m $inputMeshFile -o $OurputDir --prec $precondition_strength --log-detail 1 --dt $time_step --obj $objective --max-iter 2001 --ww 0 --wa 0
```

###### option

* `-j`:    the sub-routine  to run

* `-m`:   input surface mesh file

* `-o`:   The  ouput directory

* `--prec`:   precondition strength in the article, default to 1

* `--log-detail`:   output temporary results during the optimization.  

  > 1 : output basic results:  the asymptotic elastic tensor,  optimized surface,  objective and time step history;
  >
  > 2 : output the mesh before and after surgery 
  >
  > 3 :  output meshes in each iteration

* `--dt`:  maximal step, default to 0.1

* `--ww`:  regularization term, weight of  Willmore flow , default to 0.

* `--wa`:  regularization term, weight of  mean curvature flow , default to 0.

* `--obj`:  optimization objective,  using the following basic grammar

  > ###### Predefined item:
  >
  > * `c00`, `c01`, ...., `c55`:  components in the (asymptotic) elastic matrix (the Voigt representation of  asymptotic elastic tesnor)
  > * `bulk` :  bulk modulus
  > * `Yxyz`: Young's modulus in the direction of $(x,y,z)$ .
  > * `iso` :  isotropic penalty
  > * `S{e11,e22,e33,e23,e13,e22}`: optimize ADS under strain given in the brace. 
  >
  > ###### Combination rule:
  >
  > Linear combination of Predefined objective items is supported,  in the following format:
  >
  > ```
  > --obj c00+0.4bulk+iso+S{1,2,3,4,5,6}
  > ```

#### Sampling periodic meshes

We sample periodic meshes based on several trigonometric representations of TPMS, you call invoke it by

```bash
minsurf -j samplemesh --sample-type 12 --nsample 1 -o $outdir --sample-rand 0.5 --remesh-len 0.03 --remesh-iter 5 --remesh-smthiter 0
```

> * `--sample_type` denotes the TPMS type
> * `--nsample` denotes the number of surfaces to be sampled
> * `--sample-rand` denotes the strength of randomness (perturbation to TPMS)

#### Remesh a surface to produce periodic boundaries

Some meshes are geometrically periodic but topologically not, i.e., the vertices on opposite periodic bounaries do not overlapp through a periodic translation. 

We call the following command remeshing it to produce periodic boundaries

```bash
-j periodremesh -m $inputmesh -o $outdir --remesh-len 0.03 --remesh-smthiter 0 --remesh-iter 20
```

