import itertools

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

seq = read_fasta('rosalind_grph.txt')
snames = list(seq.keys())
tita = itertools.combinations(snames, 2)
#print(*tita, sep='\n')

def overlap_graph(nodes_list, k):
    result = []
    for vertice in nodes_list:
        if seq[vertice[0]][:k] == seq[vertice[1]][-k:]:
            result.append(vertice[::-1])
        elif seq[vertice[1]][:k] == seq[vertice[0]][-k:]:
            result.append(vertice)
        else:
            pass
    return result

result = overlap_graph(tita, 3)
for pair in result:
    print(pair[0], pair[1])
#print (len(result))
