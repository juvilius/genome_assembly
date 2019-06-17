#from graphviz import Digraph
import numpy as np

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

k = 90

class sss:
    def __init__(self, file):
        self.file = file
        sequences = set()
        seq = read_fasta(self.file)
        for sequence in seq:
            sequences.add(seq[sequence])
        self.sequences = sequences

class bruijn(sss):
    def __init__(self, file):
        sss.__init__(self, file)
        self.vertices = set()
        self.edges = []
        #f = Digraph('de_bruijn', filename='br.gv')
        #f.attr('node', shape='circle')
        for ss in self.sequences:
            for j in range(len(ss)-k):
                self.vertices.add(ss[j:j+k])
                self.vertices.add(ss[j+1:j+k+1])
                self.edges.append((ss[j:j+k], ss[j+1:j+k+1]))
        self.setedges = set(self.edges)
        print('vertices: ', len(self.vertices), '\nedges: ', len(self.edges), len(self.setedges))
        #for e in set(self.edges):
        #    f.edge(e[0], e[1])
        #f.view()


#rosa = bruijn('3_chplor.fasta')


class matrix(bruijn):
#MATRIX w/ named indices and methods
    def __init__(self, file):
        bruijn.__init__(self, file)
        self.mat = np.zeros((len(self.vertices), len(self.vertices)), dtype=np.int8)
        print('matrix initiated: ', self.mat.nbytes/1e6, 'MB', self.mat.shape)
        #print(0 == self.mat[50][50])

    def name(self):
        #na2nu name to number
        #nu2na number to name
        self.na2nu = dict()
        for na, nu in enumerate(self.vertices):
            self.na2nu[nu] = na

    def fill(self):
        for edge in self.setedges:
            self.mat[self.na2nu[edge[0]]][self.na2nu[edge[1]]] += 1
        np.set_printoptions(threshold=np.inf)
        writ = str(self.mat)
        with open('matrix_output.txt', 'w') as file:
            file.write(writ)


mat = matrix('3_chplor.fasta')
mat.name()
mat.fill()
