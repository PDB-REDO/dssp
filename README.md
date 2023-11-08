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

You can use the special branch `dssp-builder` to build dssp including all the
dependencies:


```console
git clone https://github.com/PDB-REDO/dssp.git -b dssp-builder
cd dssp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local
cmake --build build
```

That will build the mkdssp executable in your $HOME/.local directory.

But you can of course also build the classic way using the `trunk` branch in
which case you will have to have all dependencies installed first.

Usage
-----

See [manual page](doc/mkdssp.md) for more info. Or even better, see the [DSSP website](https://pdb-redo.eu/dssp).
