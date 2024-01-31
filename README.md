DSSP 4.4
========

This is a rewrite of DSSP, now offering full mmCIF support. The difference
with previous releases of DSSP is that it now writes out an annotated mmCIF
file by default, storing the secondary structure information in the
`_struct_conf` category.

Another new feature in this version of DSSP is that it now defines
Poly-Proline helices as well.

The DSSP program was designed by _Wolfgang Kabsch_ and _Chris Sander_ to
standardize secondary structure assignment. DSSP is a database of secondary
structure assignments (and much more) for all protein entries in the Protein
Data Bank (PDB). DSSP is also the program that calculates DSSP entries from
PDB entries.

DSSP does **not** predict secondary structure.

Requirements
------------

The tools are based on [libcif++](https://github.com/PDB-REDO/libcifpp)
and the code is written in modern C++ so you need a compiler capable
of handling C++17 code.

Building
--------

First make sure the following dependencies are installed:

* [zlib](https://github.com/madler/zlib) A compression library. You can also use your package manager to install this, on Debian/Ubuntu this is done using `apt-get install zlib1g-dev`
* [libmcfp](https://github.com/mhekkel/libmcfp.git) A library for parsing command line arguments
* [libcif++](https://github.com/pdb-redo/libcifpp.git) A library for reading and writing PBDx files.

When these are installed, you can then execute the following commands to build and install dssp:

```console
git clone https://github.com/pdb-redo/dssp.git
cd dssp
cmake -S . -B build
cmake --build build
cmake --install build
```

Usage
-----

See [manual page](doc/mkdssp.md) for more info. Or even better, see the [DSSP website](https://pdb-redo.eu/dssp).
