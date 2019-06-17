import random
import time
def read_fasta(file):
    seqs = ''
    for line in open(file, 'r'):
        if line.startswith('>'):
            sname = line.replace('>', '').replace('\n', '')
        else:
            sseq = line.replace('\n', '')
            seqs = seqs + sseq
    return seqs

def slice(seq, n, outfile):
    '''n is coverage'''
    index_f = 2
    leng = len(seq)
    print(leng)
    writ = '>read_0000001' + '\n' + seq[-100:] + seq[:100] + '\n'
    with open(outfile, 'a') as file:
        file.write(writ)
    for i in range(n):
        begi = 0
        while True:
            r = random.randint(300, 900)
            if (begi + r) < leng:
                writ = '>read_' + str(index_f).zfill(7) + '\n' + seq[begi:(begi+r)] + '\n'
                #print(writ, '\n', begi, r, index_f)
                begi += r + 1
                index_f += 1
                with open(outfile, 'a') as file:
                    file.write(writ)
            elif (begi + r) >= leng:
                writ = '>read_' + str(index_f).zfill(7) + '\n' + seq[begi:] + '\n'
                #print(writ, '\n', begi, r, index_f)
                index_f += 1
                with open(outfile, 'a') as file:
                    file.write(writ)
                break
            else:
                print('hm...')
                break


outfile = 'full_chplor.fasta'
oryza = read_fasta('oryza_chloroplast.fasta')
oryza = oryza
slice(oryza, 3, outfile)
#print(oryza)
