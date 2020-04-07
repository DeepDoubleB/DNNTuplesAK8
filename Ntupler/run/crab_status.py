import os, sys
#from termcolor import colored, cprint

dirs = os.listdir("crab_projects")

for d in dirs:
	os.system('crab status crab_projects/'+d)
	continue
	o = os.popen('crab status crab_projects/'+d).read().split("\n")
	for i, l in enumerate(o):
		if l.startswith("CRAB project directory:"): print l
		if l.startswith("Jobs status"):
			for j in range(5):
				if len(o[i+j]) < 2:continue
				if  any(s in o[i+j] for s in ['unsubmitted', 'idle', 'finished','running','transferred', 'transferring', 'failed']):
					print o[i+j]
			
		
