import os
from argparse import ArgumentParser

#parser = ArgumentParser()
#arser.add_argument("--include", help="Path to training sample", default=)
#args=parser.parse_args()
#print(args.include)
list_dir = 'lists_102v3'
f_dir = '102RE'

#names = ["QCD", "Glu", "Bul", "GluGluH"]
#names = ["GluGluH"]
#names = ["Hcc", "Hbb", "Z"]
#names = ["Hbb"]
names = ["multimass", 'ggHbb', 'ggHcc', "qcd", "Zqq"]

dcap = "dcap://grid-dcap-extern.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/anovak/{}".format(f_dir)
srm = "srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms/store/user/anovak/{}".format(f_dir)
fnal = False

dirs = os.popen("gfal-ls " + srm).read().split("\n")
#else:	dirs = os.popen("lcg-ls "+srm).read().split("\n")
dirs = [d.split("/")[-1] for d in dirs]

print("Processsing directories starting with: ", names)
print("=======================================")
all_files = []
print(dirs)
for i, d in enumerate(dirs):
    dir_files = []
    go = False
    for name in names:
        if d.startswith(name): go = True
    if not go: continue
    if True:
        sd = os.popen("gfal-ls " + srm + "/" + d).read().split('\n')
        #sd = os.popen("lcg-ls "+"srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms/store/user/anovak/94xv3"+"/"+d).read().split('\n')
        for s in sd:
            if not s.endswith(tuple(["madgraph", "pythia8"])): continue
            fs = []
            #print s
            sd2 = os.popen("gfal-ls " + srm + "/" + d + "/" +
                           s).read().replace("\n", "")
            sd3 = os.popen("gfal-ls " + srm + "/" + d + "/" + s + "/" +
                           sd2).read().replace("\n", "")
            sd4 = os.popen("gfal-ls " + srm + "/" + d + "/" + s + "/" + sd2 +
                           "/" + sd3).read().replace("\n", "")
            files = os.popen("gfal-ls " + srm + "/" + d + "/" + s + "/" + sd2 +
                             "/" + sd3 + '/' + sd4).read().split("\n")
            path = dcap + "/" + d + "/" + s + "/" + sd2 + "/" + sd3 + "/" + sd4
            for f in files:
                if not f.endswith("root"): continue
                fs.append(path + "/" + f)
            #print "gfal-ls "+srm+"/"+d+"/"+s+"/"+sd2+"/"+sd3+'/'+sd4
            #print files
            #s = s.replace("/pnfs/physik.rwth-aachen.de/cms", "")
            #srm = srm.replace(".de/cms/store/user/anovak/94xv3", ".de/cms")
            #sd2 = os.popen("lcg-ls "+srm+s).read().replace("\n", "").replace("/pnfs/physik.rwth-aachen.de/cms", "")
            #sd3 = os.popen("lcg-ls "+srm+sd2).read().replace("\n", "").replace("/pnfs/physik.rwth-aachen.de/cms", "")
            #sd4 = os.popen("lcg-ls "+srm+sd3).read().replace("\n", "").replace("/pnfs/physik.rwth-aachen.de/cms", "")
            #files = os.popen("lcg-ls "+srm+sd4).read().split("\n")
            #files = [f.replace("/pnfs/physik.rwth-aachen.de/cms", "") for f in files]
            #for f in files:
            #	if len(f) < 5: continue
            #	if not f.endswith("root"): continue
            #	fs.append("dcap://grid-dcap-extern.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/"+f)
            #s = s.split("/")[-1]
            #print fs

            with open(list_dir+"/" + s + ".txt", 'w') as f:
                f.write('\n'.join(fs))
