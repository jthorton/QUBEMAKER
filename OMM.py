#!/usr/bin/env python
from collections import defaultdict

import sys

# OpenMM Imports
import simtk.openmm as mm
from simtk.openmm import KcalPerKJ
import simtk.openmm.app as app 
import parmed as pmd
# ParmEd Imports
from parmed import load_file, unit as u
from parmed.charmm import CharmmParameterSet
from parmed.openmm import StateDataReporter, energy_decomposition_system
params = CharmmParameterSet('QUBE_FF_openmm.par')
ala5_gas =  load_file('new_ionized.psf')
ala5_crds = load_file('ionized.pdb')

system = ala5_gas.createSystem(params, nonbondedMethod=app.NoCutoff,nonbondedCutoff=500.0*u.angstroms, switchDistance=496.0*u.angstroms)

#serializing the system
from simtk.openmm import XmlSerializer
serialized_system = XmlSerializer.serialize(system)
outfile = open('serialized_system.xml','w')
outfile.write(serialized_system)
outfile.close()

# Create the integrator to do Langevin dynamics
integrator = mm.LangevinIntegrator(
                        300*u.kelvin,       # Temperature of heat bath
                        1.0/u.picoseconds,  # Friction coefficient
                        2.0*u.femtoseconds, # Time step
)





platform = mm.Platform.getPlatformByName('CPU')
simulation = app.Simulation(ala5_gas.topology, system, integrator, platform)
simulation.context.setPositions(ala5_crds.positions)
print ('energy from openmm library')
print (simulation.context.getState(getEnergy=True).getPotentialEnergy()) 
#platform.setPropertyValue(simulation.context, property='Precision', value='double')
# Minimize the energy
#print('Minimizing energy')

#Minimize(simulation,iters=0)

struct=pmd.load_file('ionized.pdb')
ecomps=(pmd.openmm.energy_decomposition_system(struct, system))
tot_ene=0.0
for i in range(0,len(ecomps)):
        tot_ene+=ecomps[i][1]
        print(ecomps[i][0],ecomps[i][1])
print('Total-energy %6.6f'%tot_ene)
