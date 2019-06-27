import numpy as np
import time
import random
import argparse
#from graphviz import Digraph

def argparser():
    parser = argparse.ArgumentParser(description='de Bruijn genome assembly')
    subparser = parser.add_subparsers()

    parser_c = subparser.add_parser('create', help='create artificial reads from a given sequence')
    parser_c.add_argument('-i', '--infile', required=True, type=argparse.FileType('r'), help='file with genome sequence')
    parser_c.add_argument('-n', '--coverage', type=int, default=6, help='coverage, default is 6')
    parser_c.add_argument('-l', '--length', type=int, default=100, help='variable length of reads, default is 100')
    parser_c.add_argument('-o', '--outfile', type=argparse.FileType('w'), help='fasta file with reads')

    parser_a = subparser.add_parser('assemble', help='tries to assemble reads to genome sequence')
    parser_a.add_argument('-i', '--infile', required=True, type=argparse.FileType('r'), help='file with reads')
    parser_a.add_argument('-k', '--kmer', type=int, default=43, help='k for kmers, default is 6')
    parser_a.add_argument('-p', '--origfile', required=True, type=argparse.FileType('r'), help='original sequence.fasta to proof the result')
    parser_a.add_argument('-o', '--outfile', default='results.txt', type=argparse.FileType('w'), help='result file')

    parser_d = subparser.add_parser('demo', help='takes a genome sequence, deconstructs it and tries to put it back again. for demonstration purposes.')
    parser_d.add_argument('-i', '--infile', dest='genome', required=True, type=argparse.FileType('r'), help='file with genome sequence')
    parser_d.add_argument('integer', nargs='*', type=int, action='append')
    #args = parser.parse_args(['assemble', '-i', 'p_reads_100_5.txt', '-k', '39', '-p', 'pUC19.fasta'])
    #args = parser.parse_args(['demo', '-i', 'pUC19.fasta', '70000', '7', '2000'])
    #args = parser.parse_args(['create', '-i', 'bla.txt', '-n', '5', '-l', '39', '-o', 'pUC19.fasta'])
    args = parser.parse_args()
    #print(args)
    return args


def read_fasta(file):
    seqs = {}
    for line in open(file, 'r'):
        if line.startswith('>'):
            sname = line.replace('>', '').replace('\n', '')
            seqs[sname] = ''
        else:
            sseq = line.replace('\n', '')
            seqs[sname] = seqs[sname] + sseq
    return seqs

def create_reads(infile, l, n, outfile):
    oryza = read_fasta(infile)
    seq = [oryza[sequence] for sequence in oryza][0]
    #print(seq)
    leng = len(seq)
    with open(outfile, 'w') as file:
        file.write('')
    index_f = 1
    n = int(n * leng / l)
    print('>> creating ', n, 'reads...')
    for i in range(n):
        r = random.randint(int(l - l/10), int(l + l/10))
        t = random.randint(0, leng-1)
        if leng < t + r:
            writ = '>read_' + str(index_f).zfill(7) + '\n' + seq[t:] + seq[:(r-(leng-t))] + '\n'
        else:
            writ = '>read_' + str(index_f).zfill(7) + '\n' + seq[t:(t+r)] + '\n'
        with open(outfile, 'a') as file:
            file.write(writ)
        index_f += 1

#create_reads('pUC19.fasta', 200, 5, 'p_reads_200_5.txt')




class bruijn:
    def __init__(self, file, k):
        self.file = file
        self.k = k
        sequences = set()
        seq = read_fasta(self.file)
        for sequence in seq:
            sequences.add(seq[sequence])
        self.sequences = sequences

    def kmerise(self, k):
        self.vertices = set()
        self.edges = []
#        f = Digraph('de_bruijn', filename='br.gv')
#        f.attr('node', shape='circle')
        for ss in self.sequences:
            for j in range(len(ss)-k):
                self.vertices.add(ss[j:j+k])
                self.vertices.add(ss[j+1:j+k+1])
                self.edges.append((ss[j:j+k], ss[j+1:j+k+1]))
        self.setedges = set(self.edges)
        print('vertices: ', len(self.vertices), '\nedges: ', len(self.edges), len(self.setedges))
#        for e in set(self.edges):
#            f.edge(e[0], e[1])
#        f.view()

class matrix(bruijn):
#MATRIX w/ named indices and methods
    def __init__(self, file, k):
        bruijn.__init__(self, file, k)
        self.kmerise(k)
        self.mat = np.zeros((len(self.vertices), len(self.vertices)), dtype=bool)
        print('matrix initiated: ', self.mat.nbytes/1e6, 'MB', self.mat.shape)
        self.genome = ''
        self.contigs = set()

    def name(self):
        #na2nu name to number
        #nu2na number to name
        self.na2nu = dict()
        for na, nu in enumerate(self.vertices):
            self.na2nu[nu] = na
        self.nu2na = {y:x for x,y in self.na2nu.items()}

    def fill(self):
        for edge in self.setedges:
            self.mat[self.na2nu[edge[0]]][self.na2nu[edge[1]]] = 1


    def start_path(self):
        self.nodd = np.count_nonzero(np.sum(self.mat, axis=0) - np.sum(self.mat, axis=1))
        self.nolink = len(self.vertices) - np.count_nonzero(np.sum(self.mat, axis=0) + np.sum(self.mat, axis=1))
        print('odd vertices: ', self.nodd)
        if self.nodd != 0:
            self.odd_vertices = np.sum(self.mat, axis=0) - np.sum(self.mat, axis=1)
            elementa = random.choice(np.nonzero(self.odd_vertices)[0])
            closed = 0
        else:
            elementa = random.randint(0, len(self.vertices)-1)
            closed = 0
        self.genome += self.nu2na[elementa]
        while np.count_nonzero(self.mat) > closed:
            if np.count_nonzero(self.mat[elementa]) > 0:
                #print(count, np.count_nonzero(self.mat[elementa]))
                elementb = np.nonzero(self.mat[elementa])[0][0]
                #print(elementa, elementb)
                self.mat[elementa][elementb] = 0
                self.genome += self.nu2na[elementb][-1]
                elementa = elementb
            else:
                self.contigs.add(self.genome)
                self.genome = ''
                self.nodd = np.count_nonzero(np.sum(self.mat, axis=0) - np.sum(self.mat, axis=1))
                if self.nodd != 0:
                    self.odd_vertices = np.sum(self.mat, axis=0) - np.sum(self.mat, axis=1)
                    elementa = random.choice(np.nonzero(self.odd_vertices)[0])
                    closed = 1
                else:
                    elementa = random.randint(0, len(self.vertices)-1)
                    closed = 0
                self.genome += self.nu2na[elementa]
        self.contigs.add(self.genome)

    def assemble(self):
        self.name()
        self.fill()
        self.start_path()


    def proof(self, orig_file, resultfile):
        oryza = read_fasta(orig_file)
        self.oryza = [oryza[sequence] for sequence in oryza][0]
        leng = len(self.oryza)
        self.contigswn = list()
        for cont in self.contigs:
            n = str(self.oryza + self.oryza[:self.k+3]).find(cont[:self.k+1])
            self.contigswn.append((n, cont))
        self.contigswn.sort()
        #print(*self.contigswn, sep='\n')
        # write result to file + alignment
        writ2 = ['-']*leng

        for i,s in self.contigswn:
            #print(leng, ':', len(s))
            if (i + len(s)) > leng:
                s1 = s[:(leng-i)]
                for j,c in enumerate(s1):
                    writ2[i+j] = c
                s2 = s[(leng-i):]
                for j,c in enumerate(s2):
                    writ2[j] = c
            else:
                for j,c in enumerate(s):
                    writ2[i+j] = c
        writ = ''
        for bla in writ2:
            writ += bla
        equa = 0
        for i, c in enumerate(writ):
            if c == self.oryza[i]:
                equa += 1
        ident = equa / leng
        writ = self.oryza +'\n' + writ
        with open(resultfile, 'w') as file:
            file.write(writ)
        print('identity: ', ident*100, '%')

args = argparser()
if hasattr(args, 'coverage'):
    #print(args.infile.name, args.outfile.name)
    create_reads(args.infile.name, args.length, args.coverage, args.outfile.name)
elif hasattr(args, 'kmer'):
    mat = matrix(args.infile.name, args.kmer)
    mat.assemble()
    mat.proof(args.origfile.name, args.outfile.name)
else:
    print('input error')



#mat = matrix('p_reads_100_2.txt', 43)
#mat.assemble()
#mat.proof('pUC19.fasta', 'result.txt')
