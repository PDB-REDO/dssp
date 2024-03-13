[![github CI](https://github.com/pdb-redo/dssp/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/pdb-redo/dssp/actions)
[![GitHub License](https://img.shields.io/github/license/pdb-redo/dssp)](https://github.com/pdb-redo/dssp/LICENSE)

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

A good, modern compiler is needed to build the mkdssp program since it uses
many new C++20 features.

Building
--------

The new makefile for dssp will take care of downloading and building all requirements
automatically. So in theory, building is as simple as:

```console
git clone https://github.com/PDB-REDO/dssp.git
cd dssp
cmake -S . -B build
cmake --build build
cmake --install build
```

Usage
-----

See [manual page](doc/mkdssp.md) for more info. Or even better, see the [DSSP website](https://pdb-redo.eu/dssp).
