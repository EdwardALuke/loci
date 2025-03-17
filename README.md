
# What is Loci
Loci (pronounced low-sigh) is a both a C++ library and a programming framework
specifically designed for developing computational simulations of physical
fields, such as computational fluid dynamics. One advantage the framework
provides is automatic parallelization. Once an application is described within
the Loci framework, the application can be executed in parallel without change,
even though the description within the framework included no explicit parallel
directives. A particular advantage of the programming framework is that it
provides a formal framework for the development of simulation knowledge-bases
using *logic-relational* abstractions.

The Loci framework implements declarative programming model that is
very different from traditional imperative programming languages like
C++ or Fortran.  In declarative programming models, the programmer
describes the properties of program elements and the system assembles
the components in order to satisfy some goal.  A similar set of
programming paradigms would be represented by tools such as rules in
Makefiles or SQL query languages for databases. The advantage of this
type of system is that it can often free the programmer from the
burden of planning how to coordinate computations and as a result can
eliminate many common sources of program errors.  


# Installation

## Obtaining Source Code

Start by navigating to the directory where the Loci repo will be
located (for example, $HOME/code) and clone it. After this, navigate
to the root directory and initialize the submodules:

```bash
cd $HOME/code
git clone https://github.com/EdwardALuke/loci.git
cd loci
git submodule init
```

## Dependencies
**Required:**
* MPI implementation, such as [Open MPI](https://www.open-mpi.org/) or [MPICH](https://www.mpich.org)

  This will usually be installed on most computational clusters and is
  required to access the distributed computing facilities of Loci.
  Loci should be able to work with most versions of MPI.  (If only
  interested in compiling for pre/post processing tools, it is also
  possible to compile without MPI.)
      
* [HDF5](https://www.hdfgroup.org/solutions/hdf5/)*

  Loci uses the HDF5 (Hierarchical Data Format) library for binary file
  I/O including finite volume grid formats.  Loci cannot be run
  without HDF5.

**Optional:**
* Partitioning library - Either PT-Scotch or ParMETIS
    * [PT-Scotch](https://www.labri.fr/perso/pelegrin/scotch/)*
    * ParMETIS implemented in the libraries:
      * [GKlib](https://github.com/KarypisLab/GKlib)*
      * [METIS](https://github.com/KarypisLab/METIS)*
      * [ParMETIS](https://github.com/KarypisLab/ParMETIS)*

* [PETSc](https://petsc.org/)*

  A package that allows Loci to use preconditioned Krylov subspace
  solvers for solving large sparse linear systems.
       
* [libxml2](https://github.com/GNOME/libxml2)

  A library for parsing XML formatted input files.  This library is
  used to define geometric region adaptation for the mesh
  adaptation tool.  If the library is not found, this feature
  is disabled.  Most systems will have the libxml2 library
  installed as part of the operating system.  Used for XML
  based region specification for offline mesh adaptation

* [CGNS](https://cgns.github.io/)

  Used for pre- and post-processing tools for converting from CGNS based unstructured mesh formats.


*These dependencies are available from within the Loci repository as
submodules. The latest version confirmed to work with the currently
checked-out version of Loci can be obtained via:

```bash
git submodule update ext/${modulename}
```

Follow each package's installation procedure and ensure the resulting binary
and library directories are added to your `$PATH` and `$LD_LIBRARY_PATH`
environment variables.

## Compiling Loci

Next, execute the `configure` script to setup the build directory and
files.  `configure --help` will provide the most common options. If no
build directory is provided as an argument the `OBJ` directory will
automatically be created to store the files created during
compilation. The `OBJ` directory will be the place where Loci is
temporarily installed while building Loci.  Below is a set of example
commands to configure and install Loci:

```bash
./configure --no-sub-dir --prefix=$HOME/local/
make install
```

The result of the above commands will be installed in the user-defined
*prefix* destination directory. This directory should be referenced as
`$LOCI_BASE` when linking and compiling other software built on this
Loci framework.  Note, that parallel make can be used to speedup the
compilation process using the -j <num jobs> command if you are
compiling on a multi-core system.

### Useful configure command options

The ./configure script has many useful options that can be used to
control how Loci is configured for compilation.  

* `--no-sub-dir`

  By default the configure script will install Loci into a
  sub-directory of the specified prefix directory that is given unique
  name based on the release version and compiler used.  This can be
  useful in managing multiple compiled versions on a single system.
  When the `--no-sub-dir` is specified, then Loci is configured to be
  installed directly into the prefix directory.

* `--nompi`

  This option allows compilation of a version of Loci that does not
  use the MPI library.  This can sometimes be helpful to allow an
  installed version that can run pre or post processing utilities on
  cluster head nodes where MPI execution may be disabled.  Note, this
  will require a serial version of the HDF5 libraries that are
  compiled without enabling parallel I/O features.

* `--obj-dir`

  This command takes the directory of the object directory where the
  code will be compiled.  By default this is the OBJ directory created
  in the source tree, however, it can be set to a directory such as
  /tmp/OBJ which can be useful on systems with limited storage under
  the users home directory.

* `--help`

  Provides documentation on all of the available configure options that may be useful in customizing your Loci install.


### Configuration Files

When the configure script is executed it sets up two configuration
files that adapt the Loci the specific configuration of the system.
By default, these files are placed in the OBJ directory where most of
the compilation takes place.  Most of the time the `./configure`
script will do a good job setting up the configuration, but on rare
occasions it may be necessary to edit these configuration files to
customize your Loci compile.  These files are as follows:

* `comp.conf` - (Compiler config) Used to create compiler specific flags
  (i.e. gcc, icc, etc.) that are compatible with Loci during the build process.

* `sys.conf` - (System config) Created by the `configure` script that 
  automatically finds the needed libraries and executables from the user's
  `$PATH` and `$LD_LIBRARY_PATH` environment variables or command line options.

* `Loci.conf` - A generic configuration file that sets up GNU make
  rules that are helpful in compiling Loci programs.

When Loci is installed these configuration files are copied to the
installed directory allowing application programs to be compiled with
the same setup as the Loci framework.  To take the greatest advantage
of these features, base your application makefile on the generic
Makefile provided in the `Tutorial/heat/Makefile` directory.  This
Makefile is setup such that by editing the target executable name and
specifying the `LOCI_BASE` environment variable to point to the Loci
installed directory one can easily compile Loci applications.

## For Developers

The following make commands are provided after successfully executing `./configure`:

| Command                       | Description                                                  |
| -----------                   | -----------                                                  |
| `make default`                | Compile and install Loci within the OBJ directory            |
| `make install`                | Compile and install Loci in the specified prefix directory   |
| `make uninstall`              | Undo the install (remove Loci from the install directory)    |
| `make test`                   | Run unit and regression tests                                |
| `make docs`                   | Compile latex documentation to PDF files                     |
| `make clean`                  | Clean compiled files from OBJ directory                      |
| `make distclean`              | Remove all configuration and compilation artifacts           |

## See Also

Please read the install.txt file for more documentation on Loci installation.


# Disclaimer

The source code in the sprng is the parallel random number generator library
developed at the Florida State University. It is provided here as a convenience.
The library can be found at http://sprng.cs.fsu.edu/

The source code under FVMtools/libadf is the ADF library is an open
source library that is used by some grid converters found in FVMtools
and is provided here as a convenience.

Note that the parMETIS tools are licensed only for non-commercial use by
University of Minnesota.  For commercial use, contact University of
Minnesota to obtain a license, or alternatively use the open source
PTScotch graph partitioning library.

All other source code in this directory is provided by Mississippi
State University as free software under the Lesser GNU General Public
License version 3.  See the file [](COPYING.LESSER) for license terms.

Ed Luke  (luke@cse.msstate.edu) | Professor (retired) of Computer Science and Engineering |
Mississippi State University
