# Example for: model.rename_segments()

# Assegno alle CHAINS del modello nuove lettere; rinumero i RESIDUES delle catene come erano nel PDB della XRAY.
# Assign new segment names and write out the new model:
#R= recettore, A-C= asymmetric unit of Coat Protein, G,U= Calcium ions
# rinumero in base ai PDB originari

from modeller import *

# Read the MODEL with all HETATM and water records (so there are two 'chains'):
env = environ()
env.io.atom_files_directory = ['../atom_files']
env.io.hetatm = True
env.io.water = True

mdl = model(env, file='CPabcCooP-FABP3.B99990002.pdb')
mdl.rename_segments(segment_ids=('R', 'A', 'B', 'C', 'G', 'U'),
                    renumber_residues=[1, 101, 101, 67, 1, 1])
mdl.write(file='CPCooP-FABP3.docked02.pdb')
