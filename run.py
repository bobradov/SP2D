#!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python 



# Run script

import CombinatorialLoop as cl
import SplitTable        as st
import numpy as np
import os
from subprocess import call
from parm_subs  import parm_subs

	
# Default values
default_dict = { 'L' : 5.0, 
			 'W' : 5.0, 
			 'vsp' : 2.4923, 
			 'Vg': 0.6, 
			 'Eps': 22.2,
			 'dev_type' : 'fin',
			 'lattice_thick'            : 0.5,
			 'period_thickness_ratio'   : 7.0,
			 'barrier'                  : 0.0,
			 'd2_shift'                 : 0.0
			  }
			 
# Make a combinatorial loop over a few parameters	
c_pgaa = cl.CombinatorialLoop( { 'W'   : [3.0, 5.0, 7.0],
                                 'Vg'  : np.linspace(-0.2, 0.8, 30),
                                 'd2_shift' : [0,  -0.05],
                                 'dev_type' : ['pgaa']
                                 })             
                                 
c_nw = cl.CombinatorialLoop( { 'W'        : [3.0, 4.0, 5.0, 7.0, 11.0],
                               'Vg'       : np.linspace(-0.2, 0.8, 30),
                               'dev_type' : ['nw']
                             })    
                                 
c_fin = cl.CombinatorialLoop( { 'W'        : [30.0],
                                'Vg'       : np.linspace(-0.2, 0.8, 30),
                                'dev_type' : ['fin']
                              })     
                               
split_table = st.SplitTable( [c_pgaa, c_fin], default_dict ) 
n_nodes  = 60
n_splits = 0
offset   = 0



# Create run directories       
for split, cur_split in enumerate(split_table):

	print cur_split
	
	# Make a run directiory
	split_dir = 'Split_' + str(split+offset)
	if not os.path.exists(split_dir):
		os.makedirs(split_dir)
    	
    # Put substituted version of template into run directory
	parm_subs('template.py', split_dir + '/' + 'sim.py', cur_split)
	
	n_splits += 1
	
	
	
# Create jobfiles
for cur_jobfile_index in range(0, n_nodes):
	jf_name = 'jobfile_' + str(cur_jobfile_index)
	jobf = open(jf_name, 'w')
	
	# Header
	jobf.write('# Job Name\n' ) 
	jobf.write('#BSUB -J pGAA_SP\n') 
	jobf.write('# merge stdout and stderror\n')
	jobf.write(
	'#BSUB -oo /project/SAS_xfer/data_hub/bobradovic/5nmCMOS/pGAA/pGAA_Template/OutLog.o\n')
	jobf.write(
	'#BSUB -eo /project/SAS_xfer/data_hub/bobradovic/5nmCMOS/pGAA/pGAA_Template/OutLog.e\n')
	jobf.write('#BSUB -q "advmpi.q advsmp.q"\n')
	jobf.write('#BSUB -n 16\n') 
	jobf.write('#BSUB -R "span[hosts=1]"\n\n')
	
	# Now loop over all split directories which belong to this node
	for cur_split in range(0, n_splits):
		if cur_split % n_nodes == cur_jobfile_index:
			jobf.write('cd Split_' + str(cur_split+offset) + '\n')
			jobf.write('if [ -e "FOM.dat" ]\n')
			jobf.write('then\n')
			jobf.write('echo "Split already finished." > job.stdout\n')
			jobf.write('cd ..\n')
			jobf.write('else\n')
			jobf.write(
			'/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python sim.py > job.stdout'
			+ '\n')
			jobf.write('cd ..' + '\n')
			jobf.write('fi\n')
	jobf.write('exit' + '\n')
			
	jobf.close()
	
# Now run all the jobs
for job in range(0, n_nodes):
	print 'Submitting job: ', job
	jobname = 'jobfile_' + str(job)
	os.system("bsub < " + jobname)
	


