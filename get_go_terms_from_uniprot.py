#!/usr/bin/env python

import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)

  #Input/output files
  p.add_option('-i', '--input', help='Input UniProt file [None, REQD]')
  p.add_option('-o', '--output', help='Output file with format "locus_tag<TAB>GO,GO,GO" [None,REQD]')

  opts, args = p.parse_args()

  term_dict = {}

  with open(opts.input,"r") as fin:
    name = ""
    for line in fin:
      if "OrderedLocusNames" in line:
        if name:
          term_dict[name] = terms
        names = line[5:].strip().split(";")
        for namestring in names:
          if "OrderedLocusNames" in namestring:
            name = namestring.split("=")[1].split(" ")[0]
            terms = []
      elif "GO:" in line:
        terms.append(line.strip().split(";")[1][1:])
    term_dict[name] = terms
    
  with open(opts.output, "w") as fout:
    for name in sorted(term_dict):
      fout.write("%s\t%s\n" % (name, ",".join(term_dict[name])))

if __name__=="__main__":
  main()