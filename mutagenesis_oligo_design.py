#!/usr/bin/python


from Bio.Seq import Seq
from Bio.Alphabet import generic_dna 

import optparse, os


def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-i', '--input', help='Input. Should have: fasta file, name of sequence, region start, region end. Tab-delimited [None, REQD]')
    p.add_option('-o', '--output', help='Output fasta [None, REQD]')
    p.add_option('-r', '--replacement', help='Input replacement file. Should have: Target codon, replacement codon. Tab-delimited [None,REQD]')

    opts, args = p.parse_args()

    seq_info_dict=get_seq_info(opts.input)

    oligos_unmodified=[]

    replacement_dict = {'GCC': 'GCT','GCT': 'GCC','GCA': 'GCC','GCG': 'GCC','AGA': 'AGG','AGG': 'AGA','CGG': 'AGA','CGC': 'AGA','CGA': 'AGA','CGT': 'AGA','AAC': 'AAT','AAT': 'AAC','GAC': 'GAT','GAT': 'GAC','TGC': 'TGT','TGT': 'TGC','CAG': 'CAA','CAA': 'CAG','GAG': 'GAA','GAA': 'GAG','GGC': 'GGA','GGA': 'GGC','GGG': 'GGC','GGT': 'GGC','CAC': 'CAT','CAT': 'CAC','ATC': 'ATT','ATT': 'ATC','ATA': 'ATC','CTG': 'CTC','CTC': 'CTG','TTG': 'CTG','CTT': 'CTG','TTA': 'CTG','CTA': 'CTG','AAG': 'AAA','AAA': 'AAG','ATG': 'XATGX','TTC': 'TTT','TTT': 'TTC','CCC': 'CCT','CCT': 'CCC','CCA': 'CCC','CCG': 'CCC','AGC': 'TCC','TCC': 'AGC','TCT': 'AGC','TCA': 'AGC','AGT': 'AGC','TCG': 'AGC','TGA': 'TGA','TAA': 'TAA','TAG': 'TAG','ACC': 'ACA','ACA': 'ACC','ACT': 'ACC','ACG': 'ACC','TGG': 'XTGGX','TAC': 'TAT','TAT': 'TAC','GTG': 'GTC','GTC': 'GTG','GTT': 'GTG','GTA': 'GTG'}

    translation_dict = {'GCC': 'A', 'GCT': 'A', 'GCA': 'A', 'GCG': 'A', 'AGA': 'R', 'AGG': 'R', 'CGG': 'R', 'CGC': 'R', 'CGA': 'R', 'CGT': 'R', 'AAC': 'N', 'AAT': 'N', 'GAC': 'D', 'GAT': 'D', 'TGC': 'C', 'TGT': 'C', 'CAG': 'Q', 'CAA': 'Q', 'GAG': 'E', 'GAA': 'E', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GGT': 'G', 'CAC': 'H', 'CAT': 'H', 'ATC': 'I', 'ATT': 'I', 'ATA': 'I', 'CTG': 'L', 'CTC': 'L', 'TTG': 'L', 'CTT': 'L', 'TTA': 'L', 'CTA': 'L', 'AAG': 'K', 'AAA': 'K', 'ATG': 'M', 'TTC': 'F', 'TTT': 'F', 'CCC': 'P', 'CCT': 'P', 'CCA': 'P', 'CCG': 'P', 'AGC': 'S', 'TCC': 'S', 'TCT': 'S', 'TCA': 'S', 'AGT': 'S', 'TCG': 'S', 'TGA': 'Z', 'TAA': 'Z', 'TAG': 'Z', 'ACC': 'T', 'ACA': 'T', 'ACT': 'T', 'ACG': 'T', 'TGG': 'W', 'TAC': 'Y', 'TAT': 'Y', 'GTG': 'V', 'GTC': 'V', 'GTT': 'V', 'GTA': 'V'}

    most_common_dict = {'A': 'GCC', 'R': 'AGA', 'N': 'ACC', 'D': 'GAC', 'C': 'TGC', 'Q': 'CAG', 'E': 'GAG', 'G': 'GGC', 'H': 'CAC', 'I': 'ATC', 'L': 'CTG', 'K': 'AAG', 'M': 'ATG', 'F':'TTC', 'P':'CCC', 'S':'AGC', 'Z': 'TGA', 'T':'ACC', 'W': 'TGG', 'Y':'TAC', 'V':'GTG'}

    start_coords=[]
    end_coords=[]

    # In theory this could be extended to handle multiple input sequences, but that hasn't happened yet
    for seq in seq_info_dict:
       i=int(seq_info_dict[seq][1])-15

       while i <= (int(seq_info_dict[seq][2])-15):
           oligos_unmodified.append(get_seq(seq, seq_info_dict[seq][0], i, i+32))
           start_coords.append(i)
           end_coords.append(i+32)
           i=i+3

    j=0

    for oligo in oligos_unmodified:
	target = oligo[15:18]
        following = oligo[18:21]

        pre_start = 12
        pre_end = 15
        while True:
            preceding = oligo[pre_start:pre_end]
            pre_start = pre_start-3
            pre_end = pre_end-3
            if (preceding.upper() != 'ATG' and preceding.upper() != 'TGG'):
                pre_start = pre_start+3
                pre_end = pre_end+3
                break

        preceding_modified=replacement_dict[preceding.upper()]

        foll_start=18
        foll_end=21
        while True:
            following=oligo[foll_start:foll_end]
            foll_start=foll_start+3
            foll_end=foll_end+3
            if (following.upper() != 'ATG' and following.upper() != 'TGG'):
                foll_start = foll_start-3
                foll_end = foll_end-3
                break

        following_modified=replacement_dict[following.upper()]

        target_modified = []

        for amino_acid in most_common_dict:
           if translation_dict[target.upper()] != amino_acid:
              target_modified.append(most_common_dict[amino_acid])


        for new_target in target_modified:
            new_oligo=oligo[:pre_start]+preceding_modified+oligo[pre_end:15]+new_target+oligo[18:foll_start]+following_modified+oligo[foll_end:]
            print "%s_%s_%s_%s_to_%s\t%s" % (seq_info_dict[seq][0], start_coords[j], end_coords[j], target.upper(), new_target, new_oligo.upper())

        j=j+1

def get_seq_info(bfile):
    to_change_dict={}
    fin=open(bfile, 'r')
    for line in fin:
        cols=line.strip().split('\t')
        to_change_dict[cols[0]]=[x for x in cols[1:]]
    return to_change_dict

def get_seq(file, seq_name, start, stop):
    seq_dict = read_fasta_dict(file)
    seq_of_interest = seq_dict[seq_name][int(start)-1:int(stop)]
    return seq_of_interest

def read_fasta_lists(file):
    fin = open(file, 'r')
    count=0

    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:])
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)

    return names, seqs

def read_fasta_dict(file):
    names, seqs = read_fasta_lists(file)
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict

if __name__ == '__main__':
    main()
