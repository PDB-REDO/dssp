# NAME

mkdssp - Assign secondary structure to proteins

# SYNOPSIS

mkdssp \[OPTION\] input \[output\]

# DESCRIPTION

The DSSP program was designed by Wolfgang Kabsch and Chris Sander to
standardize secondary structure assignment. DSSP is a database of
secondary structure assignments (and much more) for all protein entries
in the Protein Data Bank (PDB). mkdssp is the program that calculates
DSSP entries from PDB entries. mkdssp does **not** predict secondary
structure.

The original DSSP program wrote output in a fixed format, this version
by default writes annotated mmCIF files, storing the secondary structure
information in the \_struct_conf category.

Since version 4.0 the mkdssp program also assigns PPII helices.

# OPTIONS

The input file can be either mmCIF or PDB format and the file may be
gzip compressed. Note that input files must be formatted correctly. E.g.
PDB files must have a CRYST1 record. More info:
https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html#CRYST1

The output is optional, if omitted the output is written to *stdout*. If
the name of the output file ends with either *.gz* or *.bz2* the output
is compressed accordingly.

**\--output-format**=\[dssp\|mmcif\]

:   If an output file is specified, the extension of the filename is
    used to choose to output format, but if it is unclear, mmcif is the
    default. Use this option to force output in either the old fixed
    column DSSP format or the new annotated mmCIF format.

**\--no-dssp-categories**

:   When writing mmCIF files, suppress the output of all dssp\_
    categories.

**\--min-pp-stretch**

:   This option can be used to define the minimal number of residues
    with PHI/PSI angles within the range required to assing a PP helix.

**\--write-other**

:   By default the new format does not write the structure information
    for OTHER. Use this flag to change that.

**\--components**

:   The knowledge of compounds is loaded from the CCD file
    *components.cif* that should have been installed by *libcifpp*. You
    can override that file by using this option.

**\--extra-compounds**

:   As an addition to the standard *components.cif* file, you can add
    more files using this option. Files should be either in CCD format
    or should be CCP4 restraints files.

**\--mmcif-dictionary**

:   The default mmCIF dictionary file is installed by the *libcifpp*
    library but you can override it using this option.

# DETAILS

The DSSP algorithm assings secondary structure based on the energy
calculated for H-bonds.\
**Table 1. Secondary Structures recognized**

  --------------- -------------- -------------
     DSSP Code      mmCIF Code    Description
         H         HELX_RH_AL_P   Alphahelix
         B             STRN       Betabridge
         E             STRN         Strand
         G         HELX_RH_3T_P     Helix_3
         I         HELX_RH_PI_P     Helix_5
         P         HELX_LH_PP_P   Helix_PPII
         T          TURN_TY1_P       Turn
         S             BEND          Bend
     ' ' (space)      OTHER          Loop
  --------------- -------------- -------------

# BUGS

The mmCIF format currently lacks a lot of information that was available
in the old format like information about the bridge pairs or the span of
the various helices recognized. Also the accessibility information is
left out.

If you think this information should be part of the output, please
contact the author.

# AUTHOR

Written by Maarten L. Hekkelman \<maarten@hekkelman.com\>

# REPORTING BUGS

Report bugs at https://github.com/PDB-REDO/dssp/issues
