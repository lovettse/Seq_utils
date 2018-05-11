dat <- read.table("PBMC-K562_97664_Sept_29_2017_KA_BJH.report.csv", skip=8, header=TRUE, sep=",")

out_dat <- data.frame(row.names=unique(dat$Type))

for(cell_type in unique(dat$Type)) {
  out_dat[cell_type,"Cell_count"] <- legth(dat[dat$Type==cell_type,"NoContam.Read.Pairs"])
  out_dat[cell_type,"Reads"] <- sum(dat[dat$Type==cell_type, "NoContam.Read.Pairs"])
  out_dat[cell_type,"Reads.per.cell"] <- sum(dat[dat$Type==cell_type, "NoContam.Read.Pairs"])/length(dat[dat$Type==cell_type,"NoContam.Read.Pairs"])
  out_dat[cell_type,"Assigned"] <- sum(dat[dat$Type==cell_type, "Assigned"])
}