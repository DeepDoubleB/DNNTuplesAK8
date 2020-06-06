import os, sys
#from termcolor import colored, cprint

dirs = os.listdir("crab_projects")

for d in dirs:
	if "MX" not in d: continue
	#if "1000to1400" in d: continue
	#if "Inf" in d: continue
	#os.system('crab resubmit crab_projects/'+d+" --sitewhitelist=T2_DE_RWTH")
	#os.system('crab resubmit crab_projects/'+d)
#	continue
	print 'crab status crab_projects/'+d
	os.system('crab status crab_projects/'+d)
	#o = os.popen('crab status crab_projects/'+d).read().split("\n")
	print o
	continue
	for i, l in enumerate(o):
		if l.startswith("CRAB project directory:"): print l
		if l.startswith("Jobs status"):
			for j in range(5):
				if len(o[i+j]) < 2:continue
				if  any(s in o[i+j] for s in ['unsubmitted', 'idle', 'finished','running','transferred', 'transferring', 'failed']):
					print o[i+j]
			
		
