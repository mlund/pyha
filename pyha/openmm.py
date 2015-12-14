from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit

def findForce(system, forcetype, add=True):
  """ Finds a specific force in the system force list - added if not found."""
  for force in system.getForces():
    if isinstance(force, forcetype):
      return force
  if add==True:
    system.addForce(forcetype())
    return findForce(system, forcetype)
  return None

def setGlobalForceParameter(force, key, value):
  for i in range(force.getNumGlobalParameters()):
    if force.getGlobalParameterName(i)==key:
      print('setting force parameter', key, '=', value)
      force.setGlobalParameterDefaultValue(i, value);

def atomIndexInResidue(residue):
  """ list of atom index in residue """
  index=[]
  for a in list(residue.atoms()):
    index.append(a.index)
  return index

def getResiduePositions(residue, positions):
  """ Returns array w. atomic positions of residue """
  ndx = atomIndexInResidue(residue)
  return np.array(positions)[ndx]

def uniquePairs(index):
  """ list of unique, internal pairs """
  return list(combinations( range(index[0],index[-1]+1),2 ) )

def addHarmonicConstraint(harmonicforce, pairlist, positions, threshold, k):
  """ add harmonic bonds between pairs if distance is smaller than threshold """
  print('Constraint force constant =', k)
  for i,j in pairlist:
    distance = unit.norm( positions[i]-positions[j] )
    if distance<threshold:
      harmonicforce.addBond( i,j,
          distance.value_in_unit(unit.nanometer),
          k.value_in_unit( unit.kilojoule/unit.nanometer**2/unit.mole ))
      print("added harmonic bond between", i, j, 'with distance',distance)

def addExclusions(nonbondedforce, pairlist):
  """ add nonbonded exclusions between pairs """
  for i,j in pairlist:
    nonbondedforce.addExclusion(i,j)

def rigidifyResidue(residue, harmonicforce, positions, nonbondedforce=None,
    threshold=6.0*unit.angstrom, k=2500*unit.kilojoule/unit.nanometer**2/unit.mole):
  """ make residue rigid by adding constraints and nonbonded exclusions """
  index    = atomIndexInResidue(residue)
  pairlist = uniquePairs(index)
  addHarmonicConstraint(harmonic, pairlist, pdb.positions, threshold, k)
  if nonbondedforce is not None:
    for i,j in pairlist:
      print('added nonbonded exclusion between', i, j)
      nonbonded.addExclusion(i,j)

