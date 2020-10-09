## Bit-Exact ECC Recovery (BEER): 

This software implements BEER as described in [our MICRO 2020 academic
paper](https://people.inf.ethz.ch/omutlu/pub/BEER-bit-exact-ECC-recovery_micro20.pdf)
[1]. BEER uses a SAT solver to identify the parity-check matrix (i.e., H-matrix)
of a linear block code using only uncorrectable error patterns. BEER is useful
when dealing with an invisible or unknown ECC mechanism (e.g., modern DRAM chips
that use on-die ECC). 

[Our
paper](https://people.inf.ethz.ch/omutlu/pub/BEER-bit-exact-ECC-recovery_micro20.pdf)
[1] shows that knowing the parity-check matrix of an ECC code allows predicting
the ECC code's post-correction error characteristics, enabling:

- System designers to make informed decisions when using third-party products with unknown ECC codes.
- Test and validation engineers to diagnose the root-causes for observed post-correction errors and craft test patterns that target pre-correction error properties.
- Scientific researchers to tie post-correction error characteristics to those of the pre-correction errors.

To learn about BEER and its potential applications in much greater detail, please see any of:
- [Our paper](https://people.inf.ethz.ch/omutlu/pub/BEER-bit-exact-ECC-recovery_micro20.pdf)
- [15-minute conference talk](https://www.youtube.com/watch?v=D97oAbCaJWk)
  (slides in
  [pptx](https://people.inf.ethz.ch/omutlu/pubBEER-bit-exact-ECC-recovery_micro20-talk.pptx)
  or
  [pdf](https://people.inf.ethz.ch/omutlu/pub/BEER-bit-exact-ECC-recovery_micro20-talk.pdf))
- [90-second lightning talk](https://www.youtube.com/watch?v=hgSziiRTUY4)
  (slides in
  [pptx](https://people.inf.ethz.ch/omutlu/pub/BEER-bit-exact-ECC-recovery_micro20-lightning-talk.pptx)
  or
  [pdf](https://people.inf.ethz.ch/omutlu/pub/BEER-bit-exact-ECC-recovery_micro20-lightning-talk.pdf))

*Please send questions to Minesh Patel at minesh.patelh@gmail.com*

## BEER High-Level Overview

[Our
paper](https://people.inf.ethz.ch/omutlu/pub/BEER-bit-exact-ECC-recovery_micro20.pdf)
[1] provides a detailed overview of BEER's purpose and design, and we highly
recommend the interested user to refer to it when delving deeply into BEER.

At a high level, BEER consists of three phases:
1. BEER induces uncorrectable errors in a real DRAM chip by pausing DRAM refresh
   and allowing data-retention errors to occur (in this software, we simulate
   doing so using the [EINSim DRAM error-correction
   simulator](https://github.com/CMU-SAFARI/EINSim) [2, 3]). When inducing
   data-retention errors, BEER uses carefully-chosen test patterns that restrict
   errors to specific codeword bit positions. This is possible due to the
   inherent asymmetry of data-retention errors: data-retention errors typically
   occur only in bits whose values are programmed to specific hardware-dependent
   values (i.e., such that the underlying DRAM cell capacitor is set to its
   fully-charged state). Designing and using these test patterns is motivated
   and discussed in our paper [1].
2. For each test pattern, BEER aggregates the bit-positions of uncorrectable
   errors that were observed. The particular bit-positions depend on the
   specific parity-check matrix that the ECC code uses.
3. Using a SAT solver, BEER computes the unique parity-check matrix that can
   cause the observed uncorrectable error pattern.

This software implements the aforementioned three steps and can optionally use
the [EINSim simulator](https://github.com/CMU-SAFARI/EINSim) [2, 3] to
substitute the functionality of real experiments with physical DRAM chips.

## Codebase Overview

BEER is a C++ command-line tool set up as a Makefile project. All source files
are contained within the ```src/``` directory and library dependencies are
provided within ```lib/```.

We use Doxygen to document the source code and provide a Doxyfile for building
HTML and LaTeX documentation. To build the documentation, simply issue:

```
$ doxygen
```
when in the project directory, or point ```doxygen``` to the provided Doxyfile.
The HTML documentation will be built under ```doxygen/html/index.html```. 

## Dependencies

BEER relies on the [EINSim simulator](https://github.com/CMU-SAFARI/EINSim) [2,
3] to simulate injecting data-retention errors. EINSim must be built separately
and its executable path provided to BEER via the command line.

## Toolchain

Building and running BEER requires a working C++11 toolchain (e.g., GCC, Clang,
MSVC). BEER has been built and tested with:

	- GCC 9.1.0
	- GCC 7.4.0
	- GCC 6.3.0
	- Apple LLVM 10.0.0 (clang-1000.11.45.5)
	- Apple LLVM 9.1.0 (clang-902.0.39.2)

## Building in Linux

```
$ make [-j <# threads>] [other make options] <target>
```

The makefile has various targets, described as follows:

- ```release``` builds ```beer``` with full optimizations
- ```debug``` builds ```beer.d``` with no optimization and debug symbols
- ```all``` builds both ```release``` and ```debug```
- ```doc``` builds ```doxygen``` documentation using the provided doxyfile
- ```clean``` cleans build and binary files for both ```release``` and ```debug```

Omitting the ```target``` argument defaults to the ```release``` configuration. 

## Usage

BEER runs as a command-line tool with several CLI options that are shown when running BEER without options:

```
$ ./path/to/beer
```

BEER has two primary modes that are configured using the ```-m``` flag:
- ```-m f``` file mode, where ECC code parameters are read from a JSON configuration file provided using the ```-f``` option.
- ```-m g``` generation mode, where BEER generates a random ECC code using EINSim with the ECC code parameters provided on the CLI.

For example, the following command line will first generate a random (7,4)
Hamming code, will simulate injecting errors within the 1- and 2-CHARGED test
patterns (discussed in our paper [1]) using EINSim, and will apply BEER to
determine the original ECC function using the Z3 SAT solver:
```
./path/to/beer -m g -k 4 -p 1 -p 2 -l ALL_T path/to/einsim
```

## Licensing

The current version of the simulator is provided as-is under the MIT license.

The following header-only libraries are used and are located under ```lib``` with their own license:
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [cxxopts](https://github.com/jarro2783/cxxopts)
- [rapidjson](https://rapidjson.org/)

BEER requires the Z3 solver built as a library. ```lib``` includes the Z3-v4.8.7
source and a build script to help build Z3 as a library in-directory. However,
you may modify the Makefile to link against a different version (e.g.,
system-wide installation).
- [Z3 Solver](https://github.com/Z3Prover/z3)

## Attribution

Please cite the following paper when using BEER:

[\[1\] Minesh Patel, Jeremie S. Kim, Taha Shahroodi, Hasan Hassan, and Onur Mutlu, "Bit-Exact ECC Recovery (BEER): Determining DRAM On-Die ECC Functions by Exploiting DRAM Data Retention Characteristics", in the Proceedings of the 53rd Annual ACM/IEEE International Symposium on Microarchitecture (MICRO 2020), Virtual, October 2020.](https://people.inf.ethz.ch/omutlu/pub/BEER-bit-exact-ECC-recovery_micro20.pdf)

Other references:

[\[2\] Minesh Patel, Jeremie S. Kim, Hasan Hassan, and Onur Mutlu, "Understanding and Modeling On-Die Error Correction in Modern DRAM: An Experimental Study Using Real Devices", in the Proceedings of the 49th Annual IEEE/IFIP International Conference on Dependable Systems and Networks (DSN 2019), Portland, OR, USA, June 2019.](https://people.inf.ethz.ch/omutlu/pub/understanding-and-modeling-in-DRAM-ECC_dsn19.pdf)

[\[3\] EINSim Simulator on GitHub](https://github.com/CMU-SAFARI/EINSim)