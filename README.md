DSSP 4.3
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

You also need [libmcfp](https://github.com/mhekkel/libmcfp.git)

Building
--------

Make sure you install [libcif++](https://github.com/PDB-REDO/libcifpp) and [libmcfp](https://github.com/mhekkel/libmcfp.git) first before building.

After that, building should be as easy as typing:

```bash
git clone https://github.com/PDB-REDO/dssp.git
cd dssp
mkdir build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build
```

Usage
-----

See [manual page](doc/mkdssp.md) for more info. Or even better, see the [DSSP website](https://pdb-redo.eu/dssp).
