import os, sys
#from termcolor import colored, cprint

dirs = os.listdir("crab_projects")

for d in dirs:
    os.system('crab status crab_projects/'+d)
			
		
