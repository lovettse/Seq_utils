#!/usr/bin/env python

import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input from samogitia [None,REQD]')
  p.add_option('-o', '--output', help='Output [None,REQD]')
  
  opts, args = p.parse_args()

  out_lines = []
  
  sections = get_sections(opts.input)
    
  with open(opts.input, "r") as fin:
    header = fin.readline().strip()
    for line in fin:
      transitions_dict = zero_transition_dict(sections)
      cols = line.strip().split("\t")
      state = cols[0]
      transition_count = cols[1]
      transitions = cols[2:]

      for transition in transitions:
        transition_cols = transition[1:-1].split(",")
        source = transition_cols[2]
        dest = transition_cols[3]

        transitions_dict[(source,dest)] += 1
      new_out_line = "%s\t%s" % (state, transition_count)
      for transition_type in sorted(transitions_dict):
        new_out_line += "\t%s" % ( transitions_dict[transition_type])
      out_lines.append(new_out_line)

  transition_header = ""
  for transition in sorted(zero_transition_dict(sections)):
    transition_header += "\t%s" % ("-".join(transition))

  with open(opts.output,"w") as fout:
    fout.write("state\ttransition_count%s\n" % (transition_header))
    for line in out_lines:
      fout.write("%s\n" % line)

def zero_transition_dict(sections):
  transitions_dict = {}

  for source in sections:
    for dest in sections:
      if source != dest:
        transitions_dict[(source, dest)] = 0
  return transitions_dict

def get_sections(in_file):
  sections = set()
  
  with open(in_file,"r") as fin:
    fin.next()
    for line in fin:
      cols = line.strip().split("\t")
      info_cols = cols[2:]
      for info_col in info_cols:
        this_cols = info_col[1:-1].split(",")
        sections.add(this_cols[2])
        sections.add(this_cols[3])

  return sections

if __name__ == '__main__':
  main()