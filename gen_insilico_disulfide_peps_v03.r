# Code for the script gen_insilico_disulfide_peps_v03.r script
#
# MS/MS spectra analysis of interpeptide disulfide bridge analysis; 
# the script can take as input a list of tryptic peptides and obtain in silico
# the MS/MS tables with all theoretical b adn y ions of pairs between two
# tryptic peptides.
# 
# This script reads a txt/tab file containing a list of tryptic peptides and 
# exports the the PeptideTheoreticalIonsTable from the two peptide pairs;
# the script is based on Xu H, Zhang L, Freitas MA., J Proteome Res. 2008,
# Jan;7(1):138-44. doi: 10.1021/pr070363z.
# Requirements: 'OrgMassSpecR' package
# Plase install the package using 'install.package()' function if this is 
# not yet available in your R installation
#
#=====================gen_insilico_disulfide_peps_v03.r========================
# Prepare the working enviroment 
rm(list=ls())
# ! /Library/Frameworks/R.framework/Resources
# Enable command arguments
args=commandArgs(TRUE)
x=args[1] # directory cotaining '.tab' file of the peptides to be analyzed
# y=args[2] # .mgf file of MS/MS spectra against to search for
message('Code written by CristianV.A.Munteanu ©2025')
###############################################################################
# Developped at the Protein Protein Chemistry Facility of IBRA - RO
# (https://www.biochim.ro/facility-1/) by Dr. Cristian (V.A.) Munteanu
# Please contact cristian.munteanu@biochim.ro
# (http://cristianvamunteanu.eu) for any supplementary information.
#
# Input files
#
# .txt/tab file
#
# Version log
# Version 0.10: reads in a txt/tab file and exports the corresponding tables
# in the results_tables/ directory. The algorithm for b and y ions generation is 
# based on 'OrgMassSpecR' package. Currently the script is limited only to  
# combinations of tryptic peptides containing a single Cys residue. The script 
# is based on Xu H, Zhang L, Freitas  MA., J Proteome Res. 2008, 
# Jan;7(1):138-44.  doi: 10.1021/pr070363z.
# Requirements: input as exported by mgftoTab_v012.r script. 
# Version 0.20: added disulfide-bond specific fragment ions (see the script) 
# according to the Choi S. et. al, J Proteome Res. 2010 Jan;9(1):626-35. doi: 
# 10.1021/pr900771r. PMID: 19902913; laso extended the charge state of the 
# precursors and the fragments up to +20.
# Version 0.30: added option to import PTM on Met residues, thus automatically 
# converting a m residue into an oxidized Met when calculating the ions
###############################################################################

##########################
t1=Sys.time()
message(paste(Sys.time(),'Loading functions...'))

# Prepare the working enviroment
options(stringsAsFactors=F,help_type='text',warn=-1)

# Load the required custom functions

# Function to convert single aa. pep. to elemental or 3-letter code
ConvertPeptide_mod=function (sequence, output = "elements", IAA = TRUE) 
{
    peptideVector <- strsplit(sequence, split = "")[[1]]
    if (output == "elements") {
        FindElement <- function(residue) {
            if (residue == "A") 
                element <- c(C = 3, H = 5, N = 1, O = 1, S = 0)
            if (residue == "R") 
                element <- c(C = 6, H = 12, N = 4, O = 1, S = 0)
            if (residue == "N") 
                element <- c(C = 4, H = 6, N = 2, O = 2, S = 0)
            if (residue == "D") 
                element <- c(C = 4, H = 5, N = 1, O = 3, S = 0)
            if (residue == "E") 
                element <- c(C = 5, H = 7, N = 1, O = 3, S = 0)
            if (residue == "Q") 
                element <- c(C = 5, H = 8, N = 2, O = 2, S = 0)
            if (residue == "G") 
                element <- c(C = 2, H = 3, N = 1, O = 1, S = 0)
            if (residue == "H") 
                element <- c(C = 6, H = 7, N = 3, O = 1, S = 0)
            if (residue == "I") 
                element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)
            if (residue == "L") 
                element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)
            if (residue == "K") 
                element <- c(C = 6, H = 12, N = 2, O = 1, S = 0)
            if (residue == "M") 
                element <- c(C = 5, H = 9, N = 1, O = 1, S = 1)
            if (residue == "m") 
                element <- c(C = 5, H = 9, N = 1, O = 2, S = 1)
            if (residue == "F") 
                element <- c(C = 9, H = 9, N = 1, O = 1, S = 0)
            if (residue == "P") 
                element <- c(C = 5, H = 7, N = 1, O = 1, S = 0)
            if (residue == "S") 
                element <- c(C = 3, H = 5, N = 1, O = 2, S = 0)
            if (residue == "T") 
                element <- c(C = 4, H = 7, N = 1, O = 2, S = 0)
            if (residue == "W") 
                element <- c(C = 11, H = 10, N = 2, O = 1, S = 0)
            if (residue == "Y") 
                element <- c(C = 9, H = 9, N = 1, O = 2, S = 0)
            if (residue == "V") 
                element <- c(C = 5, H = 9, N = 1, O = 1, S = 0)
            if (residue == "C" & IAA == FALSE) 
                element <- c(C = 3, H = 5, N = 1, O = 1, S = 1)
            if (residue == "C" & IAA == TRUE) 
                element <- c(C = 5, H = 8, N = 2, O = 2, S = 1)
            return(element)
        }
        resultsVector <- c(C = 0, H = 0, N = 0, O = 0, S = 0)
        for (i in 1:length(peptideVector)) {
            resultsVector <- FindElement(peptideVector[i]) + 
                resultsVector
        }
        resultsVector <- resultsVector + c(C = 0, H = 2, N = 0, 
            O = 1, S = 0)
        return(as.list(resultsVector))
    }
    if (output == "3letter") {
        FindCode <- function(residue) {
            if (residue == "A") 
                let <- "Ala"
            if (residue == "R") 
                let <- "Arg"
            if (residue == "N") 
                let <- "Asn"
            if (residue == "D") 
                let <- "Asp"
            if (residue == "C") 
                let <- "Cys"
            if (residue == "E") 
                let <- "Glu"
            if (residue == "Q") 
                let <- "Gln"
            if (residue == "G") 
                let <- "Gly"
            if (residue == "H") 
                let <- "His"
            if (residue == "I") 
                let <- "Ile"
            if (residue == "L") 
                let <- "Leu"
            if (residue == "K") 
                let <- "Lys"
            if (residue == "M") 
                let <- "Met"
            if (residue == "F") 
                let <- "Phe"
            if (residue == "P") 
                let <- "Pro"
            if (residue == "S") 
                let <- "Ser"
            if (residue == "T") 
                let <- "Thr"
            if (residue == "W") 
                let <- "Trp"
            if (residue == "Y") 
                let <- "Tyr"
            if (residue == "V") 
                let <- "Val"
            return(let)
        }
        codes <- sapply(peptideVector, FindCode)
        return(paste(codes, collapse = ""))
    }
}

# Functions to obtain in silico disulfide peptides MS/MS spectra
fragmentDisulfide=function(peptideA,peptideB,fragments="by",IAA=F,N15=F,
    custom = list())
{
    results_list <- vector("list")
    for (sequence_number in 1:length(peptideA)) {
        peptide_vector <- strsplit(peptideA[sequence_number],
            split = "")[[1]]
        peptide_length <- length(peptide_vector)
        if (peptide_length < 2)
            stop("sequence must contain two or more residues")
        C <- 12
        H <- 1.0078250321
        O <- 15.9949146221
        S <- 31.97207069
        N <- ifelse(N15 == TRUE, 15.0001088984, 14.0030740052)
        proton <- 1.007276466
        electron <- 0.00054857990943
        residueMass <- function(residue) {
            if (residue == "A")
                mass = C * 3 + H * 5 + N + O
            if (residue == "R")
                mass = C * 6 + H * 12 + N * 4 + O
            if (residue == "N")
                mass = C * 4 + H * 6 + N * 2 + O * 2
            if (residue == "D")
                mass = C * 4 + H * 5 + N + O * 3
            if (residue == "E")
                mass = C * 5 + H * 7 + N + O * 3
            if (residue == "Q")
                mass = C * 5 + H * 8 + N * 2 + O * 2
            if (residue == "G")
                mass = C * 2 + H * 3 + N + O
            if (residue == "H")
                mass = C * 6 + H * 7 + N * 3 + O
            if (residue == "I")
                mass = C * 6 + H * 11 + N + O
            if (residue == "L")
                mass = C * 6 + H * 11 + N + O
            if (residue == "K")
                mass = C * 6 + H * 12 + N * 2 + O
            if (residue == "M")
                mass = C * 5 + H * 9 + N + O + S
            if (residue == "m")
                mass = C * 5 + H * 9 + N + O * 2 + S
            if (residue == "F")
                mass = C * 9 + H * 9 + N + O
            if (residue == "P")
                mass = C * 5 + H * 7 + N + O
            if (residue == "S")
                mass = C * 3 + H * 5 + N + O * 2
            if (residue == "T")
                mass = C * 4 + H * 7 + N + O * 2
            if (residue == "W")
                mass = C * 11 + H * 10 + N * 2 + O
            if (residue == "Y")
                mass = C * 9 + H * 9 + N + O * 2
            if (residue == "V")
                mass = C * 5 + H * 9 + N + O
	    if (residue == "U")
                mass = C * 3 + H * 5 + N + O + S+MonoisotopicMass(
                ConvertPeptide_mod(peptideB,IAA=F))-2*1.0078250321
            if (residue == "C" & IAA == FALSE)
                mass = C * 3 + H * 5 + N + O + S
            if (residue == "C" & IAA == TRUE)
                mass <- ifelse(N15 == FALSE, C * 5 + H * 8 +
                  N * 2 + O * 2 + S, C * 5 + H * 8 + N + 14.0030740052 +
                  O * 2 + S)
            if (length(custom) != 0)
                for (i in 1:length(custom$code)) if (residue ==
                  custom$code[i])
                  mass = custom$mass[i]
            return(mass)
        }
        masses <- sapply(peptide_vector, residueMass)
        pm <- sum(masses)
        p1 <- round(pm + H * 2 + O + proton, digits = 3)
        p2 <- round((pm + H * 2 + O + (2 * proton))/2, digits = 3)
        p3 <- round((pm + H * 2 + O + (3 * proton))/3, digits = 3)
        if (fragments == "by") {
            b1 <- vector(mode = "numeric", length = 0)
            b2 <- vector(mode = "numeric", length = 0)
            bs <- vector(mode = "character", length = 0)
            bi <- vector(mode = "integer", length = 0)
            y1 <- vector(mode = "numeric", length = 0)
            y2 <- vector(mode = "numeric", length = 0)
            ys <- vector(mode = "character", length = 0)
            yi <- vector(mode = "integer", length = 0)
            for (i in 1:(peptide_length - 1)) {
                mass <- sum(masses[1:i])
                b1[i] <- round(mass + proton, digits = 3)
                b2[i] <- round((b1[i] + proton)/2, digits = 3)
                bs[i] <- paste(peptide_vector[1:i], collapse = "")
                bi[i] <- i
            }
            for (j in 2:peptide_length) {
                mass <- sum(masses[j:peptide_length])
                y1[j - 1] <- round(mass + H * 2 + O + proton,
                  digits = 3)
                y2[j - 1] <- round((y1[j - 1] + proton)/2, digits = 3)
                ys[j - 1] <- paste(peptide_vector[j:peptide_length],
                  collapse = "")
                yi[j - 1] <- peptide_length - j + 1
            }
            ms1seq <- rep(peptideA[sequence_number], times = ((2 *
                (length(bi))) + (2 * (length(yi)))))
            ms1z1 <- rep(p1, times = ((2 * (length(bi))) + (2 *
                (length(yi)))))
            ms1z2 <- rep(p2, times = ((2 * (length(bi))) + (2 *
                (length(yi)))))
            ms1z3 <- rep(p3, times = ((2 * (length(bi))) + (2 *
                (length(yi)))))
            ms2seq <- c(rep(bs, times = 2), rep(ys, times = 2))
            b1.type <- paste("[b", bi, "]1+", sep = "")
            b2.type <- paste("[b", bi, "]2+", sep = "")
            y1.type <- paste("[y", yi, "]1+", sep = "")
            y2.type <- paste("[y", yi, "]2+", sep = "")
            ms2type <- c(b1.type, b2.type, y1.type, y2.type)
            ms2mz <- c(b1, b2, y1, y2)
        }
        if (fragments == "cz") {
            c1 <- vector(mode = "numeric", length = 0)
            c2 <- vector(mode = "numeric", length = 0)
            cs <- vector(mode = "character", length = 0)
            ci <- vector(mode = "integer", length = 0)
            z1 <- vector(mode = "numeric", length = 0)
            z2 <- vector(mode = "numeric", length = 0)
            zs <- vector(mode = "character", length = 0)
            zi <- vector(mode = "integer", length = 0)
            for (i in 1:(peptide_length - 1)) {
                mass <- sum(masses[1:i])
                c1[i] <- round(mass + 3 * H + N + proton, digits = 3)
                c2[i] <- round((c1[i] + proton)/2, digits = 3)
                cs[i] <- paste(peptide_vector[1:i], collapse = "")
                ci[i] <- i
            }
            for (j in 2:peptide_length) {
                mass <- sum(masses[j:peptide_length])
                z1[j - 1] <- round(mass + O - N, digits = 3)
                z2[j - 1] <- round((z1[j - 1] + proton)/2, digits = 3)
                zs[j - 1] <- paste(peptide_vector[j:peptide_length],
                  collapse = "")
                zi[j - 1] <- peptide_length - j + 1
            }
            ms1seq <- rep(peptideA[sequence_number], times = ((2 *
                (length(ci))) + (2 * (length(zi)))))
            ms1z1 <- rep(p1, times = ((2 * (length(ci))) + (2 *
                (length(zi)))))
            ms1z2 <- rep(p2, times = ((2 * (length(ci))) + (2 *
                (length(zi)))))
            ms1z3 <- rep(p3, times = ((2 * (length(ci))) + (2 *
                (length(zi)))))
            ms2seq <- c(rep(cs, times = 2), rep(zs, times = 2))
            c1.type <- paste("[c", ci, "]1+", sep = "")
            c2.type <- paste("[c", ci, "]2+", sep = "")
            z1.type <- paste("[z", zi, "]1+", sep = "")
            z2.type <- paste("[z", zi, "]2+", sep = "")
            ms2type <- c(c1.type, c2.type, z1.type, z2.type)
            ms2mz <- c(c1, c2, z1, z2)
        }
        results_list[[sequence_number]] <- data.frame(ms1seq,
            ms1z1, ms1z2, ms1z3, ms2seq, ms2type, ms2mz)
    }
    return(as.data.frame(do.call("rbind", results_list)))
} # modified version of FragmentPeptide function from OrgMassSpecR 

disulfidePrIonTable=function(peptideA,peptideB,fragments='by'){
 frag_table=rbind(fragmentDisulfide(peptideA=gsub('C','U',peptideA) 
 ,peptideB=peptideB,IAA=F,fragments=fragments),fragmentDisulfide( 
 peptideA=gsub('C','U',peptideB),peptideB=peptideA,IAA=F,fragments=fragments))
 frag_table[,1]=gsub('U','C',frag_table[,1])
 # Add peptide chain ID for each b/y ion series
 frag_table[,6]=apply(frag_table[c(1,6)],1,function(x) {
  if(x[1]==peptideA){
   return(paste0('A',x[2]))
  }else{
   return(paste0('B',x[2]))
  }
 })
 return(frag_table)
} # function to obtain a merged product ion table from both peptide chains

# Load the required packages
require('OrgMassSpecR') # for peptide fragmentation
##########################
message(paste(Sys.time(),'Reading input file...'))
# Read input files
peps=read.delim(x)

filename=sub('\\.[^.]*.$','',x) # get the name of the file

# Keep only unique pep. seqs.
peps=unique(peps)

message(paste(Sys.time(),'File contains',nrow(peps),'unique peptides...'))

# Assemble the table with unique combinations of peptide pairs
peps=setNames(data.frame(t(combn(peps$PepSeq,2))),c('peptideA','peptideB'))

message(paste(Sys.time(),'Obtained',nrow(peps),'pairs of peptides...'))

message(paste(Sys.time(),'Assembling product ion tables for export...'))
suppressWarnings(dir.create(paste0(filename,'_insilico_disulfide_peps_res/'))) 
# create directory to retain all tables
th_path=paste0(filename,'_insilico_disulfide_peps_res/')

# Calculate disulfide bridge peptide pairs product ions table
peps_msms=invisible(lapply(1:nrow(peps),function(i){
 # Get pep. position
 cat('\rFinished',i,'of',nrow(peps)) # export MS/MS table

 # Calculate 'b-/y-' ions m/z
 temp=data.frame(apply(disulfidePrIonTable(peptideA=peps[i,1],peptideB=peps[i
 ,2],fragments='by'),2,function(x) gsub('U','C',x))) # get the ion table
 temp[,c(2:4,7)]=apply(temp[,c(2:4,7)],2,function(x) as.numeric(x)) # conv.

 # Add 'b0-/y0-' ions m/z to the table (water loss)
 water=round(2*1.007276466+15.9949146221,3) # calc. water Mw for 'STED' aa.
 water_nl=temp # get 'b-/y-' m/z values
 water_nl[,7]=water_nl[,7]-water/as.numeric(gsub('.*]|\\+','',water_nl[,6]))
 stop_pos=(nchar(unique(temp[,1])[1])-1)*4 # get chain A eding pos.
 water_nl[,5]=c(gsub('C',unique(water_nl[,1])[2],water_nl[1:stop_pos,5])
 ,gsub('C',unique(water_nl[,1])[1],water_nl[(stop_pos+1):nrow(water_nl),5]))
 water_nl[,7]=as.numeric(apply(water_nl[,c(5,7)],1,function(x)
 if(length(grep('S|T|E|D',x[1]))==0){NA}else{x[2]})) # keep sel. 'b0/y0' ions
 water_nl[,5]=temp[,5] # adjust seqs.
 water_nl[,6]=paste0(substr(water_nl[,6],1,3-1),0,substr(water_nl[,6],3
 ,nchar(water_nl[,6])),sep='') # adjust ion ids.
 water_nl=water_nl[which(complete.cases(water_nl)),] # keep only rows of int.
 # temp=rbind(temp,water_nl) # bind the two dfs.

 # Add 'b*-/y*-' ions m/z to the table (ammonia loss)
 ammonia=round(3*1.007276466+14.0030740052,3) # calc. ammonia Mw
 ammonia_nl=temp # get 'b-/y-' m/z values
 ammonia_nl[,7]=ammonia_nl[,7]-ammonia/as.numeric(gsub('.*]|\\+',''
 ,ammonia_nl[,6])) # substract ammonia loss from 'b-/y-' ions m/z
 ammonia_nl[,5]=c(gsub('C',unique(ammonia_nl[,1])[2],ammonia_nl[1:stop_pos,5])
 ,gsub('C',unique(ammonia_nl[,1])[1],ammonia_nl[(stop_pos+1):nrow(ammonia_nl) 
 ,5])) # replace 'C' with the pep. cross-linked sequence
 ammonia_nl[,7]=as.numeric(apply(ammonia_nl[,c(5,7)],1,function(x)
 if(length(grep('R|K|N|Q',x[1]))==0){NA}else{x[2]})) # keep sel. 'b*/y*' ions
 ammonia_nl[,5]=temp[,5] # adjust seqs.
 ammonia_nl[,6]=paste0(substr(ammonia_nl[,6],1,3-1),'*',substr(ammonia_nl[,6] 
 ,3,nchar(ammonia_nl[,6])),sep='') # adjust ion ids.
 ammonia_nl=ammonia_nl[which(complete.cases(ammonia_nl)),] # keep sel. rows
 temp=rbind(temp,water_nl,ammonia_nl) # bind dfs.

 # Add 'b3+/y3+' m/z values to the product ion table
 three_charge=temp # get product ion table
 three_charge[,7]=as.numeric(apply(three_charge[,6:7],1,function(x)
 if(as.numeric(gsub('.*]|\\+','',x[1]))==1){x[2]}else{NA})) # conv. to NA m. ch.
 three_charge=three_charge[which(complete.cases(three_charge)),] # keep sing. ch
 three_charge[,7]=(three_charge[,7]+2*1.007276466)/3 # calc. +3 m/z values
 three_charge[,6]=gsub(']1',']3',three_charge[,6])
 temp=rbind(temp,three_charge) # bind the two dfs.
 rm(ammonia_nl,water_nl,stop_pos,three_charge) # clean the enviroment

 # Save the table for +1, +2 & +3 charges
 peps_msms=temp

 # Calculate m/z for fragment ions up to the +20 charge
 temp=data.frame('ms1seq'=rep(peps_msms[grep('\\]1+',peps_msms$ms2type) 
 ,'ms1seq'],length(4:20))
 ,'ms1z1'=rep(peps_msms[grep('\\]1+',peps_msms$ms2type),'ms1z1'],length(4:20)) 
 ,'ms1z2'=rep(peps_msms[grep('\\]1+',peps_msms$ms2type),'ms1z2'],length(4:20)) 
 ,'ms1z3'=rep(peps_msms[grep('\\]1+',peps_msms$ms2type),'ms1z3'],length(4:20))
 ,'ms2seq'=rep(peps_msms[grep('\\]1+',peps_msms$ms2type),'ms2seq'],length(4:20))
 ,'ms2type'=unlist(lapply(4:20,function(i){ # get the ion ids.
  gsub('\\]1+',paste0('\\]',i),peps_msms[grep('\\]1+',peps_msms$ms2type) 
  ,'ms2type'])
 }))
 ,'ms2mz'=unlist(lapply(4:20,function(i){
  return((peps_msms[grep('\\]1+',peps_msms$ms2type),'ms2mz'] 
  +1.007276466*(i-1))/i)
 })))
  
 # Merge the lower & high-charge state dfs.
 peps_msms=rbind(peps_msms,temp)
 
 # Add the m/z ms1 values for the higher charge states
 ms1_mz_temp=unlist(lapply(4:20,function(i){
  (unique(peps_msms$ms1z1)+(i-1)*1.007276466)/i
 })) # get m/z values
 names(ms1_mz_temp)=unlist(lapply(4:20,function(i){
  paste0('ms1z',i)
 })) # add the names for each charge state
 ms1_mz=data.frame(do.call(cbind,lapply(1:length(ms1_mz_temp),function(i){
  rep(ms1_mz_temp[i],nrow(peps_msms))
 }))) # convert to df.
 names(ms1_mz)=names(ms1_mz_temp) # name the columns
 rm(ms1_mz_temp) # clean the enviroment
 row.names(ms1_mz)=NULL # remove rownames annotation
 
 # Add the ms1 m/z values
 peps_msms=data.frame(cbind(peps_msms[,c('ms1seq','ms1z1','ms1z2','ms1z3')] 
 ,ms1_mz,peps_msms[,c('ms2seq','ms2type','ms2mz')]))
 
 # Export the product ion table
  # write.table(temp,paste0(filename,'_insilico_disulfide_peps_res/' 
  # ,paste0(paste0(c(peps[i,1],peps[i,2]),collapse='_'),'.tab')),sep='\t' 
  # ,row.names=F,col.names=T)
  write.table(peps_msms,substring(paste0(filename,'_insilico_disulfide_peps_res/' 
  ,paste0(paste0(c(peps[i,1],peps[i,2]),collapse='_'),'.tab')),1,175),sep='\t' 
  ,row.names=F,col.names=T)
 return(temp)
}))

message()
message(paste(Sys.time(),"Finished")) 
t2=Sys.time() 
time=t2-t1 
units(time)='secs' 
message(paste0('File processed in ',round(time,2),' secs.'))
message()
