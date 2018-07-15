#!/usr/bin/env python

import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input fasta  (must be aligned) [None,REQD]')
  p.add_option('-r', '--ref', help='Fasta entry to use as reference (must be full header) [None, REQD]')
  
  opts, args = p.parse_args()

  in_dict = read_fasta_dict(opts.input)
  
  ref = in_dict[opts.ref]
  
  for header in in_dict:
    if not header == opts.ref:
      print header
      query = in_dict[header]
      for x in range(len(in_dict[header])):
        if not query[x] == ref[x]:
          print "%s\t%s\t%s" % (ref[x], x+1, query[x])
      print ""   
          

def read_fasta_dict(file):
  names, seqs = read_fasta_lists(file)
  fasta_dict = dict(zip(names, seqs))
  return fasta_dict

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
        seqs.append(seq.strip())
      seq=''
    else:
      seq += line
  seqs.append(seq)

  return names, seqs


if __name__ == '__main__':
  main()