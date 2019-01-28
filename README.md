# QUBEMAKER

This repo contains the work in progress to implement the QUBE FF into OPENMM.

The current version of the QUBE.xml is based on the AMBER FF meaning all residue templates are the same and should be constructed in an amber consistent manner.

## Details 

### Atom Types
The QUBE FF needs every atom type to be different as we have more bond and angle terms than other FF’s I have changed all of the atom names to new atom types. This means that some bond, angle and torsion terms may be defined multiple times.


### Bonds, Angles and Torsions
Bonds are currently grouped together in residues and listed in the same order, this makes finding them easy just search for the residue tag of the term you are changing.
Angles and torsions also follow this pattern. 

When making Improper torsion entries the first atom of the dihedral MUST BE THE CENTRAL ATOM, followed by the other three atoms in alphabetical order. 

### Nonbonded
The nonbonded section is taken straight from amber just to assign some parameters to the atoms. There is then going to be a script entered into the FF xml file which will contain the new nonbonded information and when the FF is loaded the script will run and update the values to the QUBE values and also apply the opls geometric combination rules.

### QUBEMAKER
This is a python module which takes QUBEMAKER input files? and pdb and produces the specific QUBE.xml FF file. The top half of every FF is the same and the script only edits the nonbonded section. 

### Ideas
We can include QUBEMAKER into QUBEKit or make it its own python package?
 
We can make QUBEMAKER all in python as the basic QUBE.xml means we no longer need the matlab script. We can instead make something that will take the results directly from ONETEP and produce the specific xml for the protein?

I think we have to prepare the proteins using amber as the xml file is based on the AMBER topologies. 

### TODO

The rest of the bonded infomation needs filling in for every residue and the energys comparing to using the par and psf files.

Once we have a new input parameter format for the nonbonded scetion we need a script to add this to the qube.xml in the QUBEMAKER.py script.

