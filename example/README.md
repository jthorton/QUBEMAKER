# QUBEMAKER

## Example

This script contains an ideal example of how to use the QUBE FF. Here we can load the qube.xml file for the system as a normal xml file. 
Here we have to use a specific OPLS style tip3p water xml file during the solvation command due to the combination rules being different to the standard water models.
Once solvated we can write a new pdb to make sure it has been completed.
Next we have to edit the LJ combination rules for every atom in the system so we use a slightly different fix to normal.

Note in this example the charges and LJ terms to be changed are at the top of the script instead of being taken from the psf and par files. 
This will be  automated in future but first we should decided where the information will be coming from.  

Running this file should take the test.pdb file apply the qube.xml to it solvated the system using tip3p_opls.xml , minimize and then run the simulation for 10000 steps.

 