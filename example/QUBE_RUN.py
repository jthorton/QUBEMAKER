
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import parmed as pmd
import sys


class OpenMM:
    """A class to load and run a system in openmm using QUBEMAKER"""

    def __init__(self, pdb_file, xml, opls=True):
        self.pdb = pdb_file
        self.xml = xml
        self.opls = opls
        self.simulation = None
        self.system = None
        self.modeller = None
        self.natoms = 0

        self.charge = [-0.5264, 0.1504, 0.1504, 0.1504, 0.7857, -0.6082, -0.6104, 0.3534, 0.1797, 0.0778, -0.2488, 0.0848, 0.0637, 0.3273, -0.0217, -0.3691, 0.0868, 0.0868, 0.0868, -0.3739, 0.0939, 0.0939, 0.0939, 0.4661, -0.62, -0.4412, 0.3083, -0.0041, 0.0613, 0.0613, 0.0613]
        self.eps = [0.061947, 0.035964, 0.035964, 0.035964, 0.061947, 0.100097, 0.126338, 0, 0.061947, 0.035964, 0.061947, 0.035964, 0.035964, 0.061947, 0.035964, 0.061947, 0.035964, 0.035964, 0.035964, 0.061947, 0.035964, 0.035964, 0.035964, 0.061947, 0.100097, 0.139874, 0, 0.061947, 0.035964, 0.035964, 0.035964]
        self.eps = [eps*4.184 for eps in self.eps]
        self.sigma = [2.025268, 1.297126, 1.297126, 1.297126, 1.675650, 1.712584, 1.825603, 0, 1.838169, 1.282129, 1.962368, 1.380017, 1.366169, 1.724631, 1.407049, 2.004426, 1.373238, 1.373238, 1.373238, 1.991692, 1.353407, 1.353407, 1.353407, 1.841623, 1.761646, 1.763415, 0, 1.847901, 1.374892, 1.374892, 1.374892]
        print(len(self.sigma))
        print(len(self.eps))
        assert len(self.sigma) == len(self.eps)
        self.sigma = [(sigma*2**(5/6))/10 for sigma in self.sigma]

        self.openmm_system()

    def openmm_system(self):
        """Initialise the OpenMM system we will use to evaluate the energies."""

        # load the initial coords into the system and initialise
        pdb = app.PDBFile(self.pdb)

        # count the amount of atoms in the structure
        for atom in pdb.topology.atoms():
            self.natoms += 1

        # must use a custom version of the tip3p water moddel
        forcefield = app.ForceField(self.xml, 'tip3p_opls.xml')
        self.modeller = app.Modeller(pdb.topology, pdb.positions)  # set the intial positions from the pdb

        # Now we need to solvate the system
        self.modeller.addSolvent(forcefield, model='tip3p', padding=1 * unit.nanometer)

        # write out the solvated system coords
        app.PDBFile.writeFile(self.modeller.topology, self.modeller.positions, open('output.pdb', 'w'))

        # now we create the system and add the lj correction
        self.system = forcefield.createSystem(self.modeller.topology, nonbondedMethod=app.PME, constraints=None, nonbondedCutoff=1.1 * unit.nanometer)
        if self.opls:
            self.opls_lj()

        # set control parameters
        temperature = 298.15 * unit.kelvin
        integrator = mm.LangevinIntegrator(temperature, 5 / unit.picoseconds, 0.001 * unit.picoseconds)

        # add preasure to the system
        self.system.addForce(mm.MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin))

        # create the simulation context
        self.simulation = app.Simulation(self.modeller.topology, self.system, integrator)

        # set the positions to the solvated system
        self.simulation.context.setPositions(self.modeller.positions)

        # get the energy of the sytem
        print(self.simulation.context.getState(getEnergy=True).getPotentialEnergy())

        # check the energy break down
        struct = pmd.load_file('output.pdb')
        ecomps = (pmd.openmm.energy_decomposition_system(struct, self.system))
        tot_ene = 0.0
        for i in range(0, len(ecomps)):
            tot_ene += ecomps[i][1]
            print(ecomps[i][0], ecomps[i][1])
        print('Total-energy %6.6f' % tot_ene)

        # serialize the system
        serialized_system = mm.XmlSerializer.serialize(self.system)
        outfile = open('test.xml', 'w')
        outfile.write(serialized_system)
        outfile.close()

        # now minimize the system
        self.simulation.minimizeEnergy(maxIterations=100)

        # now run the simulation
        self.simulation.reporters.append(app.PDBReporter('run.pdb', 1000))
        self.simulation.reporters.append(app.StateDataReporter('run.txt', 1000, step=True, potentialEnergy=True, temperature=True))
        self.simulation.step(100000)

    def opls_lj(self, excep_pairs=None, normal_pairs=None):
        """This function changes the standard OpenMM combination rules to use OPLS, execp and normal pairs are only
        required if their are virtual sites in the molecule."""

        from numpy import sqrt

        # get system information from the openmm system
        forces = {self.system.getForce(index).__class__.__name__: self.system.getForce(index) for index in
                  range(self.system.getNumForces())}
        # use the nondonded_force tp get the same rules
        nonbonded_force = forces['NonbondedForce']
        lorentz = mm.CustomNonbondedForce(
            'epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)*4.0')
        lorentz.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        lorentz.addPerParticleParameter('sigma')
        lorentz.addPerParticleParameter('epsilon')
        lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
        self.system.addForce(lorentz)
        LJset = {}
        # Now for each particle calculate the combination list again
        for index in range(nonbonded_force.getNumParticles()):
            charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
            # print(nonbonded_force.getParticleParameters(index))
            if index < self.natoms:
                LJset[index] = (self.sigma[index], self.eps[index], self.charge[index])
                lorentz.addParticle([self.sigma[index], self.eps[index]])
                nonbonded_force.setParticleParameters(
                    index, self.charge[index], self.sigma[index], self.eps[index] * 0)
            else:
                LJset[index] = (sigma, epsilon, charge)
                lorentz.addParticle([sigma, epsilon])
                nonbonded_force.setParticleParameters(index, charge, sigma, epsilon * 0)
        for i in range(nonbonded_force.getNumExceptions()):
            (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
            # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
            # FORCE
            lorentz.addExclusion(p1, p2)
            if eps._value != 0.0:
                sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
                eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
                charge = 0.5 * (LJset[p1][2] * LJset[p2][2])
                nonbonded_force.setExceptionParameters(i, p1, p2, charge, sig14, eps)
        # for i in range(lorentz.getNumExclusions()):
        #     print(lorentz.getExclusionParticles(i))
            # If there is a virtual site in the molecule we have to change the exceptions and pairs lists
            # old method that needs updating
            # if excep_pairs:
            #     for x in range(len(excep_pairs)):  # scale 14 interactions
            #         if p1 == excep_pairs[x, 0] and p2 == excep_pairs[x, 1] or p2 == excep_pairs[x, 0] and p1 == excep_pairs[
            #             x, 1]:
            #             charge1, sigma1, epsilon1 = nonbonded_force.getParticleParameters(p1)
            #             charge2, sigma2, epsilon2 = nonbonded_force.getParticleParameters(p2)
            #             q = charge1 * charge2 * 0.5
            #             # print('charge %s'%q)
            #             sig14 = sqrt(sigma1 * sigma2) * 0.5
            #             eps = sqrt(epsilon1 * epsilon2) * 0.5
            #             nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
            # if normal_pairs:
            #     for x in range(len(normal_pairs)):
            #         if p1 == normal_pairs[x, 0] and p2 == normal_pairs[x, 1] or p2 == normal_pairs[x, 0] and p1 == \
            #                 normal_pairs[
            #                     x, 1]:
            #             charge1, sigma1, epsilon1 = nonbonded_force.getParticleParameters(p1)
            #             charge2, sigma2, epsilon2 = nonbonded_force.getParticleParameters(p2)
            #             q = charge1 * charge2
            #             # print(q)
            #             sig14 = sqrt(sigma1 * sigma2)
            #             eps = sqrt(epsilon1 * epsilon2)
            #             nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)

        return self.system



run = OpenMM('test.pdb', 'qube.xml')