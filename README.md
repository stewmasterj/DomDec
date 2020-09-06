# domdec - Domain Decomposition Preprocessor

This is an example of two types of domain decomposition methods for better load
 balancing than the simple cartesian topology suitable for distributed memory
 parallel computing of large and variable point density geometric models e.g. 
 for molecular dynamics (MD) or Peridynamics (PD).

Includes concepts of Load Balancing and Communication minimization.

Rectangular subdomains for easy communication topology.

Uses:

-Orthogonal Recursive Bisection (ORB)   
 Best for load balancing, more complicated communication setup than simple cartesian.
 Number of subdomains should be a power of 2 for optimal load balancing.

-Orthoginal Recursive Multisection (ORM)   
 Multiple paritioning, in this instance specify number of subdomains in X, Y and Z.
 Not limited to powers of 2, just `Nx*Ny*Nz = Ndomains`.

------------------

## Examples

Example of ORB using 8 domains:
![alt text](https://github.com/stewmasterj/DomDec/blob/master/screenshots/tri8ORB.png "N=8 ORB")

Example of ORM using 8 domains (4,2,1):
![alt text](https://github.com/stewmasterj/DomDec/blob/master/screenshots/tri4-2-1ORM.png "N=8 ORM")

Example of ORM using 8 domains (2,4,1):
![alt text](https://github.com/stewmasterj/DomDec/blob/master/screenshots/tri2-4-1ORM.png "N=8 ORM")

Example of ORB on a tube geometry using 8-64 domains:
![alt text](https://github.com/stewmasterj/DomDec/blob/master/screenshots/mont.png "N=4-32 ORB")


## source files

-`domdec.f90`   
 domain decomposition preprocessor, requires `ArgsMod.f90` from adjacent project.

-`testddCAF.f90`   
  a Fortran Co-Array program for testing `domdec`  

-`run.sh`  
  a shell script for running tests.

------------------

## Usage for each program

### DomDec

    Domain Decomposition for non-local particle interactions
    Usage: domdec [OPTIONS] FILE
    Options:
      -n i   ORB: Number of images or subdomains to create. [2] Default
      -N i3  ORM: Number of images along each direction. [2]
      -r r   Non-local cutoff radius. [0.867]
      -p     Use Periodic Boundary conditions. [F]
      -v     Verbose output data to stdout. [F]
      -f s   Use coarse grid file s to decompose. [FILE]
      -c i   Number of data columns in FILE. [4]
      -b r3  Override model BoxSize or supply if not in FILE
   
    Output files include:
     dd.top        topology file, data for each subdomain.
     dd.com        particle IDs to send and get from each subdomain.
     dd-out.xyz    config showing what subdomain each particle is in.
     FILE_XXXX.dd final output files to be read by simulation program.
                    this is not an XYZ file, just column data.
                    the metadata is in dd.com
   
    Input file format (XYZ)
      Number_of_points
      LX LY LZ [optional]
      XX x y z data
      ...
 The first two lines are for the number of points and boxsize, 
  then is a list of the particles, such as element name e.g. Cu
  with it's (x,y,z) coordinate and any other data about it, such
  charge or coordination, etc.

 Notes:
  1. ORB = Orthogonal Recursive Bisection
  2. ORM = Orthogonal Recursive Multisection
  3. for ORB (using option -n) best to use powers of 2 (e.g.: 1,2,4,8,16..)
  4. if you can't use a power of 2, try using ORM, PRODUCT(i(3))=NumOImages
  5. September 6th, 2020 by Ross J. Stewart

### testddCAF

    Test program for the domdec program. This uses fortran coarrays.
    Usage: cafrun -np N testddCAF rc pbc FILE
      N    is the total number of images to run this program with.
      rc   is the cutoff radius for the particle interactions, used in domdec
            to create the proper ghost domains for communication to neighbouring
            images.
      pbc  Set to T for use of periodic boundary conditions, domdec must have -p
      FILE input file basename to read with _XXX.dd extension for each image.
            for single image, FILE is just the file.
   
    Outputs
     this program produces a set of files from each image: dmp_XXX.conf
      it contains the elemnt, position, energy, flag; for each particle,
      where flag is 1 if it is a local particle and 0 if a ghost particle.

------------------

## versions

Thread model: posix
gcc version 6.3.0 20170516 (Debian 6.3.0-18+deb9u1) 

OpenCoarrays Coarray Fortran Compiler Wrapper (caf version 1.7.4)
Copyright (C) 2015-2016 Sourcery, Inc.
