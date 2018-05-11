#!/usr/local/bin/python2.7

# Sean Lovett - 10/1/15 
# Tested with SnpEff 4.1
# Modified to use database created from GFF 4/12/16

from __future__ import division
import vcf, optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]\nOutput is to stdout"
  p = optparse.OptionParser(usage)

  #Input/output files
  p.add_option('-i', '--input', help='Input VCF from SNPEff [None,required]')
  p.add_option('-m', '--min_freq',type="float", default=0, help='Minimum percent frequency to report minor variant [0]')
  p.add_option('-d', '--min_reads', type="int", default=0, help='Minimum read depth to make call [0]')
  p.add_option('-s', '--single', action="store_true", default=False, help='Set -s flag when calling variants in only one sample')

  opts, args = p.parse_args()
  
  vcf_reader = vcf.Reader(open(opts.input,'r'))

  header = ["CHROM","POS","REF","ALT","TYPE", "SEVERITY", "LOCUS", "NUC_CHANGE", "AA_CHANGE"] + vcf_reader.samples
  print "\t".join(header)

  for record in vcf_reader:
    alts=",".join(str(a) for a in record.ALT)
    info=[]
    genotypes=[]
    i=0
    if not 'ANN' in record.INFO:
      continue 

    annotation = record.INFO['ANN'][0].split("|")
    
    info.append(annotation[1])
    info.append(annotation[2])
    info.append(annotation[3])
    info.append(annotation[9][2:])
    info.append(annotation[10][2:])

    for sample in record.samples:
      if sample.data.AD: 
        total_reads = sum(sample.data.AD)

      if sample.data.GT == "." or total_reads < opts.min_reads:
        genotypes.append("--")
      else:
        out_str = []
        perc = (sample.data.AD[0]/total_reads) * 100
     

        if perc >= 100-opts.min_freq:
          out_str.append("%s" % (record.REF))
        elif perc >= opts.min_freq:
          out_str.append("%s (%.2f%%)" % (record.REF, perc))

        j=0

        for count in sample.data.AD[1:]:
          perc = int(count)/total_reads*100

          if perc >= opts.min_freq:
            out_str.append("%s (%.2f%%)" % (record.ALT[j], perc))
          j+=1

        if len(out_str) == 1 and out_str[0] == record.REF:
          genotypes.append(".")
        elif len(out_str) == 1:
          genotypes.append(out_str[0].split(" ")[0])
        else:
          genotypes.append(", ".join(out_str))

    if opts.single:
      info_str = "\t".join(str(l) for l in info)
      genotype_str = "\t".join(str(g) for g in genotypes)
      print "%s\t%d\t%s\t%s\t%s\t%s" % (record.CHROM, record.POS, record.REF, alts, info_str,genotype_str)
    else:
#      if map(str,genotypes).count(str(genotypes[0])) != len(genotypes) and map(str,genotypes).count("--") + map(str,genotypes).count(".") != len(genotypes):
      if map(str,genotypes).count("--") + map(str,genotypes).count(".") != len(genotypes):
        info_str = "\t".join(str(l) for l in info)
        genotype_str = "\t".join(str(g) for g in genotypes)
        print "%s\t%d\t%s\t%s\t%s\t%s" % (record.CHROM, record.POS, record.REF, alts, info_str,genotype_str)
  
if __name__=="__main__":
  main()
