# Homology modeling by the automodel class
# CA -- 2020
# Model: CP of TBSV plus combination of BBB crossing and MH peptides
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

# log output: verbose and all levels
log.verbose()
log.level(1, 1, 1, 1)

# Create a new MODELLER environment to build this model in
env = environ()
env.io.atom_files_directory = './:../atom_files'
env.io.hetatm = True

# Build 30 models, and assess with DOPE
a = automodel(env,
              alnfile  = 'align.ali',
              knowns   = ('CP_CooP_mod25_renum','FABP3_CPCOOP'),
              sequence = 'CPabcCooP-FABP3',
	      assess_methods=(assess.DOPE, assess.GA341)
	      )
a.starting_model= 1
a.ending_model  = 10

#a.final_malign3d= True	# all the model are fit to templates:_fit.pdb
a.deviation     = 4.0	# has to >0 if more than 1 model


# Very through VTFM optimization:
a.library_schedule = autosched.slow
a.max_var_iterations = 300

# Thorough MD optimization:
a.md_level = refine.very_slow

# Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E7
a.repeat_optimization = 2
a.max_molpdf = 1e6

a.make()                                     # do the actual homology modeling

# Ranks the models by molpdf score
ok_models = filter(lambda x: x['failure'] is None, a.outputs)
key = 'molpdf'
ok_models.sort(lambda a,b: cmp(a[key], b[key]))
m = ok_models[0]
print "Top model: %s (molpdf %.3f)" % (m['name'], m[key])
