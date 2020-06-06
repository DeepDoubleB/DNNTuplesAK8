import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-a", '--train', help="Number of files for training dataset", type=int)
parser.add_argument("-b", '--test', help="Number of files for testing dataset", type=int)
args = parser.parse_args()

# dirs = ['Hcc', 'Hbb', "qcd"]
dirs = ['qcd']
#dirs = ['multi']
maind = '/net/scratch_cms3a/novak/102REv3/'
for d in dirs:
    flist = os.listdir(maind + d)
    flist = [f for f in flist if f.endswith(".root")]
    ntot = len(flist)
    assert ntot > args.train + args.test, "Not enough files for specified split"

    f = open("list_{}.txt".format(d), "w")
    for line in flist[:args.train]:
        f.write(maind + d + "/" + line + "\n")
    f.close()

    f = open("list_test_{}.txt".format(d), "w")
    for line in flist[args.train : args.train+args.test]:
        f.write(maind + d + "/" + line + "\n")
    f.close()

