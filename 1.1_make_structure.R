# Script to make structure file

library(data.table)
library(dplyr)

# get genotypes

system("mkdir data/out/structure/temp")
system("cut -d ' ' -f 7- data/out/plink/gl_sec_het_rhina.ped > data/out/structure/temp/A")

# convert values
system("sed 's/A/1/g' data/out/structure/temp/A > data/out/structure/temp/B")
system("sed 's/T/2/g' data/out/structure/temp/B > data/out/structure/temp/C")
system("sed 's/G/3/g' data/out/structure/temp/C > data/out/structure/temp/D")
system("sed 's/C/4/g' data/out/structure/temp/D > data/out/structure/temp/E")
system("sed 's/0/-9/g' data/out/structure/temp/E > data/out/structure/temp/F")

# get ids
system("awk '{print $2}' data/out/plink/gl_sec_het_rhina.ped > data/out/structure/temp/ids")

# make structure input file

gen <- fread("data/out/structure/temp/F", header = F)

meta <- fread("data/meta/rhina_metadata.csv") %>%
  dplyr::select(c(id, pop)) %>%
  mutate(pop = case_when(grepl("oman", pop) ~ "oman",
                              grepl("saudi", pop) ~ "saudi_arabia",
                              grepl("sharjah", pop) ~ "uae",
                              grepl("abudhabi", pop) ~ "uae",
                              grepl("rak", pop) ~ "uae",
                              grepl("dubai", pop) ~ "uae",
                              grepl("srilanka", pop) ~ "srilanka",
                              grepl("bangladesh", pop) ~ "bangladesh"))

ids <- fread("data/out/structure/temp/ids", header = F) %>%
  left_join(meta, by = c("V1" = "id")) %>%
  mutate(pop = gsub(" ", "_", pop)) %>%
  mutate(POP = case_when(pop == "oman" ~ 1,
                         pop == "saudi_arabia" ~ 2,
                         pop == "uae" ~ 3,
                         pop == "srilanka" ~ 4,
                         pop == "bangladesh" ~ 5)) %>%
  mutate(popinfo = 1) %>%
  mutate(LocData = POP) %>%
  select(V1, POP, popinfo, LocData)

write.table(ids, "data/out/structure/temp/ids_pop.txt", col.names = F, row.names = F, quote = F)
system("paste -d ' ' data/out/structure/temp/ids_pop.txt data/out/structure/temp/F > data/out/structure/temp/rhina_handmade_rad.stru")

# sort structure files by population

system("sort -k 2 data/out/structure/temp/rhina_handmade_rad.stru > data/out/structure/rhina_handmade_rad_sort.stru")
system("rm -r data/out/structure/temp")

scp data/out/structure/rhina_handmade_rad_sort.stru ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/structure/data/raw/rhina_handmade_rad_sort.stru
  
# nrep <- 10
# up_to_k <- 8
# niter <- 100000
# burnin <- 50000

# ParallelStructure::parallel_structure(structure_path=STR_path,
#                                       joblist="/exports/cmvm/eddie/eb/groups/ogden_grp/structure/jobs.txt",
#                                       n_cpu=14,
#                                       infile=infile,
#                                       outpath=outpath,
#                                       numinds = nrow(gen),
#                                       numloci=(ncol(gen)-2)/2,
#                                       locdata = 0,
#                                       popdata = 1,
#                                       popflag = 0,
#                                       usepopinfo = 0,
#                                       locprior=0,
#                                       printqhat=1,
#                                       plot_output=0,
#                                       onerowperind=1)