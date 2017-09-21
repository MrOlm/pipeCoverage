#!/usr/bin/env python3

import argparse
import sys
import os
import pandas as pd

# Version 0.5
# 9/18/17
# Fixed minor error message bug
# Increase speed of gen_genome_coverage_table
# Add breadth to gen_genome_coverage_table

# Version 0.4
# 5/1/17
# Added some methods with .stb and genomes

# Version 0.3
# 1/12/17
# Added breadth functionality

# Version 0.2
# 10/15/16


class Neuron:
    def __init__(self, name, range, coverage, id, window):
        self.name = name
        self.range = range
        self.coverage = coverage
        self.id = id
        self.window = window

    def __str__(self):
        return "name: {0}, range: {1}, id: {2}".format(self.name,self.range,self.id)

    def getScaffold(self):
        words = self.name.split('_')
        return ('_'.join(words[:len(words)-1]))

class Scaffold(object):
    def __init__(self, name, length, coverage):
        self.name = name
        self.length = length
        self.coverage = coverage

class Bcov:
    def __init__(self, name, neurons, scaffolds):
        self.name = name
        self.neurons = neurons
        self.scaffolds = scaffolds

    def __str__(self):
        return "name: {0}, scaffolds: {1}, neurons: {2}".format(self.name,len(self.scaffolds),len(self.neurons))

    #return neuron names, sorted by ID
    def neuronNames(self):
        return [i for i in sorted(self.neurons,key=lambda x: self.getID(x))]

    def getNeurons(self):
        return [self.neurons[x] for x in self.neurons]

    def getID(self,neuron):
        return self.neurons[neuron].id

    def scaffoldNames(self):
        return [i for i in self.scaffolds]

    def ncov(self,neuron):
        return self.neurons[neuron].coverage

    def scov(self,scaff):
        if scaff not in self.scaffolds:
            print("{0} is not in the file {1}".format(scaffold,self.name))
        return self.scaffolds[scaff].coverage

    def getReadLength(self):
        if hasattr(self,'RL'):
            return self.RL

        if not hasattr(self,'path'):
            print("Must have path attribute to get read length")
            return None

        else:
            self.RL = parse_RL(self.path)
            return self.RL

    def getCoverageTable(self,rc=False):
        import pandas as pd

        if rc:
            RL = self.getReadLength()

        Table = {'scaffold':[],'coverage':[],'length':[]}
        if rc: Table['read_count'] = []
        for i in self.scaffolds:
            scaff = self.scaffolds[i]
            Table['scaffold'].append(scaff.name)
            Table['coverage'].append(scaff.coverage)
            Table['length'].append(scaff.length)
            if rc: Table['read_count'].append(int(round((scaff.coverage * scaff.length)/RL,0)))

        return pd.DataFrame(Table)

    def getBreadthTable(self,rc=False,min_cov=1):
        import pandas as pd

        # Calculate the breadth for each scaffold
        s2c = gen_scaff2covs(self.getNeurons())
        s2b = {s:calc_breadth(s2c[s], min_cov=min_cov) for s in s2c.keys()}

        if rc:
            RL = self.getReadLength()
        Table = {'scaffold':[],'coverage':[],'length':[],'breadth':[]}
        if rc: Table['read_count'] = []
        for i in self.scaffolds:
            scaff = self.scaffolds[i]
            Table['scaffold'].append(scaff.name)
            Table['coverage'].append(scaff.coverage)
            Table['length'].append(scaff.length)
            Table['breadth'].append(s2b[scaff.name])

            if rc: Table['read_count'].append(int(round((scaff.coverage * scaff.length)/RL,0)))

        return pd.DataFrame(Table)

    def getSplitTable(self,rc=False):
        import pandas as pd

        if rc:
            RL = self.getReadLength()
        Table = {'scaffold':[],'coverage':[],'length':[],'split':[]}
        if rc: Table['read_count'] = []

        for n in self.getNeurons():
            scaff = n.getScaffold()

            Table['scaffold'].append(scaff)
            Table['split'].append(n.name)
            Table['coverage'].append(n.coverage)
            Table['length'].append(len(n.range))

            if rc: Table['read_count'].append(int(round((n.coverage * len(n.range))/RL,0)))

        return pd.DataFrame(Table)


# group .pcov files based on what's before the "-vs-"
def group_pcovs(lis):
    groups = {}
    for path in lis:
        name = os.path.basename(path)
        base = name.split('-vs-')[0]
        groups.setdefault(base,[]).append(path)
    return groups

def parse_RL(file):
    for line in open(file).readlines():
        if line.startswith('# read length:'):
            return float(line.split()[3])

def fixname(name, id):
    if id == 'scaffold':
        scaffold_words = name.split('_')
        if 'scaffold' in scaffold_words:
            i = scaffold_words.index('scaffold')
            name = '_'.join(scaffold_words[:i+2])
        return name

    if id == 'neuron':
        words = name.split('_')
        if 'scaffold' in words:
            w = words[-1]
            i = words.index('scaffold')
            name = '_'.join(words[:i+2])
            name = name + "_" + w
        return name

def parse_bcov(cov,min_l=3000,fix=False):
    neurons = {}
    scaffolds = {}
    l = 1 # To keep track of the length range of windows
    w = 1 # To keep track of the window number
    i = 1 # To keep track of the total number of neurons
    for line in open(cov).readlines():
        if line.startswith('#'):
            continue # Header
        name,coverage,length = line.strip().split('\t')
        length = int(length)
        coverage = float(coverage)

        if not name.startswith('>'): # Scaffold
            if fix: name = fixname(name,'scaffold')
            scaffolds[name] = Scaffold(name,length,coverage)
            l = 1
            w = 1

        if name.startswith('>'): # Neuron
            if fix: name = fixname(name,'neuron')
            r = range(l,l+length)
            l += length
            if length < min_l: continue
            neurons[name[1:]] = Neuron(name[1:],r,coverage,i,w)
            i += 1
            w += 1

    B = Bcov(os.path.basename(cov),neurons,scaffolds)
    B.path = cov
    return B

def parse_bcovs(covs, min_l=0, fix=False):
    bcovs = []
    for cov in covs:
        bcovs.append(parse_bcov(cov,min_l,fix))
    return(bcovs)

def gen_scaff2covs(neurons):
    s2c = {}
    for n in neurons:
        scaff = n.getScaffold()
        s2c.setdefault(scaff,[]).append(n.coverage)
    return s2c

def calc_breadth(cov_list, min_cov=1):
    total = len(cov_list)
    b = len([x for x in cov_list if x > min_cov])
    return b/total

# Determine if all Bcovs passed in are mapping to the same assembly
def same_assembly(Bcovs):
    n_names = set(Bcovs[0].neuronNames())
    for Bcov in Bcovs:
        if n_names != set(Bcov.neuronNames()):
            print("{0} is not mapping to the same assembly as {1}".format(Bcovs[0].name,Bcov.name))
            return False
    return True

# Auto-group mappings based on the first 12 characters
def autogroup(Bcovs,out):
    groups = {}
    for B in Bcovs:
        name = out + B.name[0:12]
        if name not in groups: groups[name] = []
        groups[name].append(B)
    return(groups)

def write_output(coverage,names,learn,Bcovs,out,header):
    if coverage: write_coverage(Bcovs,out,header)
    if names: write_esomNames(Bcovs[0],out)
    if learn: write_learn(Bcovs,out)

def write_coverage(Bcovs,out,header=False,neurons=False):
    if not neurons:
        scaffolds = Bcovs[0].scaffoldNames()
        with open(out + ".cov", 'w') as o:
            if header:
                o.write("# scaffold\t")
                for B in Bcovs:
                    o.write(B.name + "\t")
                o.write("\n")
            for scaff in scaffolds:
                o.write(scaff)
                for B in Bcovs:
                    o.write("\t" + str(B.scov(scaff)))
                o.write('\n')

    else:
        neurons = Bcovs[0].neuronNames()
        with open(out + ".cov", 'w') as o:
            if header:
                o.write("# neuron\t")
                for B in Bcovs:
                    o.write(B.name + "\t")
                o.write("\n")
            for n in neurons:
                o.write(n)
                for B in Bcovs:
                    o.write("\t" + str(B.ncov(n)))
                o.write('\n')

def write_esomNames(Bcov,out):
    neurons = Bcov.getNeurons()
    neurons = sorted(neurons, key = lambda x: x.id)

    o = open(out + "_esom.names",'w')
    o.write('% {0}\n'.format(len(neurons)))
    for n in neurons:
        r = n.range
        start = r[0]
        end = r[-1]
        o.write('{0}\t{1}\t{2}:({3},{4}),\n'.format(n.id,n.name,n.getScaffold(),\
        start,end))

def write_learn(Bcovs,out):
    neurons = Bcovs[0].neuronNames()

    o = open(out + "_esom.lrn",'w')

    # Write header
    o.write("% {0}\n".format(len(neurons)))
    o.write("% {0}\n".format(len(Bcovs)+1))
    o.write("% 9\t{0}\n".format('\t'.join(['1']*(len(Bcovs)))))
    o.write("% key")
    for B in Bcovs: o.write("\t{0}".format(B.name))
    o.write("\n")

    # Write body
    for n in neurons:
        o.write(str(Bcovs[0].getID(n)) + "\t")
        for B in Bcovs:
            o.write("{0}\t".format(B.ncov(n)))
        o.write("\n")


    o.close()

def gen_genome_coverage_table(pcovs, stb, min_c = 1):
    # Load scaffold coverage info
    db = pcovs_to_df(pcovs)

    # Load scaffold to bin information
    stb = b2s_to_s2b(parse_stb(stb))
    db['bin'] = db['scaffold'].map(stb)

    # Convert to genome to bin
    Table = {'length':[],'read_count':[],'sample':[],'genome':[],'RL':[],'breadth':[]}
    for sample, x in db.groupby('sample'):
        for binn, d in x.groupby('bin'):
            if len(d) == 0:
                print('error: ',sample,bin)
                continue
            length = d['length'].sum()

            Table['length'].append(length)
            Table['read_count'].append(d['read_count'].sum())
            Table['sample'].append(sample)
            Table['RL'].append(d['RL'].tolist()[0])
            Table['genome'].append(binn)
            Table['breadth'].append(d['length'][d['coverage'] > min_c].sum() / length)

    Gdb = pd.DataFrame(Table)
    Gdb['coverage'] = (Gdb['RL'] * Gdb['read_count'])/Gdb['length']

    return(Gdb)

def pcovs_to_df(group):
    db = pd.DataFrame()
    for pcov in group:
        name = pcov.split('-vs-')[1].replace('.pcov','')
        Pcov = parse_bcov(pcov,fix=True)
        d = Pcov.getCoverageTable(rc=True)
        d['sample'] = name
        d['RL'] = Pcov.getReadLength()
        db = pd.concat([db,d])
    return(db)

def parse_stb(stb):
    bins = {}
    with open(stb, "r") as ins:
        for line in ins:
            linewords = line.strip().split('\t')
            scaffold,b = linewords[:2]
            if b not in bins:
                bins[b] = []
            bins[b].append(scaffold.strip())
    return bins

def b2s_to_s2b(b2s):
    s2b = {}
    for b in b2s:
        for s in b2s[b]:
            s2b[s] = b
    return(s2b)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    #Input Arguments
    InArgs = parser.add_argument_group('INPUT/OUTPUT')
    InArgs.add_argument('-h', action="help", help="show this help message and exit")
    InArgs.add_argument("-p", "--bcovs", nargs='*', help=".bcov file(s) (from pipe_coverage.awk)", required = True)
    InArgs.add_argument("-o", "--out", help="output basename")

    #Operational Arguments
    OpArgs = parser.add_argument_group('OPPERATIONS')
    OpArgs.add_argument("-a", "--all", help="generate all possible outputs", action = 'store_true')
    OpArgs.add_argument("-c", "--coverage", help="generate coverage file of complete scaffolds (CONCOCT format)", action = 'store_true')
    OpArgs.add_argument("-n", "--names", help="generate esom.names files", action = 'store_true')
    OpArgs.add_argument("-l", "--learn", help="generate single UN-NORMALIZED esom.lrn file", action = 'store_true')
    OpArgs.add_argument("-auto", "--auto", help="group .pcov filenames based on the first 12 characters", action = 'store_true')

    #Other Arguments
    OtherArgs = parser.add_argument_group('OTHER')
    OtherArgs.add_argument("--min_window", help="minimum window size to allow for .esom files", default = 3000)
    OtherArgs.add_argument("--header", help="write a header in coverage table (will not work natively with CONCOCT)", action='store_true')
    OtherArgs.add_argument("--fix_names", help="attempt to fix scaffold names messed up by SNAP", action='store_true')
    args = parser.parse_args()

    ### Validate input
    if args.all:
        args.coverage = True
        args.names = True
        args.learn = True
    if args.out == None:
        out = os.path.abspath('./') + '/'
    else:
        out = os.path.abspath(args.out)
        if os.path.isdir(out): out = out + '/'

    ### Parse .bcov files
    bcovs = parse_bcovs(args.bcovs,args.min_window,args.fix_names)

    ### If auto, group them into groups
    if args.auto: groups = autogroup(bcovs,out)
    else: groups = {out:bcovs}

    ### Validate groups of .bcovs that map to the same assembly
    for group in groups:
        if not same_assembly(groups[group]):
           print("\nERROR: All the input files are not mapping to the same assembly.")
           sys.exit()

    ### Write outputs
    for group in groups:
        write_output(args.coverage,args.names,args.learn,groups[group],group,args.header)
