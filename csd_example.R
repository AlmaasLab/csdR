#!/usr/bin/env Rscript
#CSD code

source("find_rho_and_var.R")

#IMPORTERE FIL OG TRANSPONERE MATRISE
# `log_progress` er importert fra find_rho_and_var.R og skriver meldingen sammen med tidspunktet
log_progress("Starting to read files")
x_1 <- read.table("healthy.txt", header = TRUE, row.names = 1) %>% as.matrix() %>% t()
x_2 <- read.table("sick.txt", header = TRUE , row.names = 1) %>%as.matrix()%>%t()
log_progress("Files imported and matrices transposed")

set.seed(123) #Sikrer reproduserbarhet
#Kjøring av CSD
csd_df <- run_csd(x_1,x_2,n_it=2000L,nThreads=5L,verbose=TRUE)

#Filter result
n_pairs <- nrow(csd_df)
pairs_to_pick <- 1000
# Vi kunne for dette formålet også godt ha skrevet
# 1:n_pairs (som jo er mer lesbart),
# men det under er en mer vanntett syntaks
# for spesielle grensetilfeller som når pairs_to_pick < n_pairs eller når n_pairs == 0
index_vector <- seq_len(min(pairs_to_pick,n_pairs))
# `order` er en indirekte sorteringsfunksjon og gir ut indeksene av den opprinnelige vektoren i den sorterte.
# For eksempel er x[order(x)] identisk med sort(x)
# Vi plukker ut de indeksene som gir høyest verdi av cVal
# Merk også at det er mer lesbart å skrive csd_df$cVal i stedet for csd_df[,7]
log_progress("Sorting C-values")
c_filter <- order(csd_df$cVal,decreasing = TRUE)[index_vector]
c_frame <- csd_df[c_filter,]

log_progress("Sorting S-values")
s_filter <- order(csd_df$sVal,decreasing = TRUE)[index_vector]
s_frame <- csd_df[s_filter,]

log_progress("Sorting D-values")
d_filter <- order(csd_df$dVal,decreasing = TRUE)[index_vector]
d_frame <- csd_df[d_filter,]

#CSD filter
# c_filter, s_filter og d_filter er nå heltallsvektorer. For å finne indeksene som inngår i alle tre gjør vi
# som følger
cs_filter <- union(c_filter, s_filter)
csd_filter <- union(cs_filter, d_filter)
csd_frame <- csd_df[csd_filter,]

#Skrive til fil 
log_progress("Writing to file")
write.table(x = c_frame, file = "C_links.txt", sep = '\t', row.names = FALSE, quote = FALSE)
write.table(x = s_frame, file = "S_links.txt", sep = '\t', row.names = FALSE, quote = FALSE)
write.table(x = d_frame, file = "D_links.txt", sep = '\t', row.names = FALSE, quote = FALSE)
write.table(x = csd_frame, file = "CSD_links.txt", sep = '\t', row.names = FALSE, quote = FALSE)
log_progress("DONE")
