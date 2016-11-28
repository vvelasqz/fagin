# subset(x, evalue > 1e-10)
# get query id column
# remove duplicates (unique)


#for P.troglodytes
Hsapiens_Blast_Ptroglodytes <- read.table("Hsapiens_Blast_Ptroglodytes")
Hsapiens_Blast_Ptroglodytes_filtered <- subset(Hsapiens_Blast_Ptroglodytes, Hsapiens_Blast_Ptroglodytes$V3 > 1e-010)
Hsapiens_Blast_Ptroglodytes_filtered_unique <- unique(Hsapiens_Blast_Ptroglodytes_filtered$V1)
write(paste0(Hsapiens_Blast_Ptroglodytes_filtered_unique), 'Hsapiens_Blast_Ptroglodytes_filtered_evalue-10_unique.txt')

#for G.gorilla
Hsapiens_Blast_Ggorilla <- read.table("Hsapiens_Blast_Ggorilla")
Hsapiens_Blast_Ggorilla_filtered <- subset(Hsapiens_Blast_Ggorilla, Hsapiens_Blast_Gorilla$V3 > 1e-010)
Hsapiens_Blast_Ggorilla_filtered_unique <- unique(Hsapiens_Blast_Ggorilla_filtered$V1)
write(paste0(Hsapiens_Blast_Ggorilla_filtered_unique), 'Hsapiens_Blast_Ggorilla_filtered_evalue-10_unique.txt')

#for N.leucogenys
Hsapiens_Blast_Nleucogenys <- read.table("Hsapiens_Blast_Nleucogenys")
Hsapiens_Blast_Nleucogenys_filtered <- subset(Hsapiens_Blast_Nleucogenys, Hsapiens_Blast_Nleucogenys$V3 > 1e-010)
Hsapiens_Blast_Nleucogenys_filtered_unique <- unique(Hsapiens_Blast_Nleucogenys_filtered$V1)
write(paste0(Hsapiens_Blast_Nleucogenys_filtered_unique), 'Hsapiens_Blast_Nleucogenys_filtered_evalue-10_unique.txt')

#for P.abelii
Hsapiens_Blast_Pabelii <- read.table("Hsapiens_Blast_Pabelii")
Hsapiens_Blast_Pabelii_filtered <- subset(Hsapiens_Blast_Pabelii, Hsapiens_Blast_Pabelii$V3 > 1e-010)
Hsapiens_Blast_Pabelii_filtered_unique <- unique(Hsapiens_Blast_Pabelii_filtered$V1)
write(paste0(Hsapiens_Blast_Pabelii_filtered_unique), 'Hsapiens_Blast_Pabelii_filtered_evalue-10_unique.txt')

#for M.mulatta
Hsapiens_Blast_Mmulatta <- read.table("Hsapiens_Blast_Mmulatta")
Hsapiens_Blast_Mmulatta_filtered <- subset(Hsapiens_Blast_Mmulatta, Hsapiens_Blast_Mmulatta$V3 > 1e-010)
Hsapiens_Blast_Mmulatta_filtered_unique <- unique(Hsapiens_Blast_Mmulatta_filtered$V1)
write(paste0(Hsapiens_Blast_Mmulatta_filtered_unique), 'Hsapiens_Blast_Mmulatta_filtered_evalue-10_unique.txt')
