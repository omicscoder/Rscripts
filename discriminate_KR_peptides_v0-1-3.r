#

# pcSILAC data analysis of 
# A375-ST cells labeled with R6/R10,K0/K0 and chased for various 
# times (t0-t5) in R0/R0,K4/K8 after EDEM2 induction with TET
# The data analysis workflow will closelly follow the one described 
# Fierro-Monti paper (Fierro-Monti et al.,PLoS One. 2013 Nov 27;8(11):e80423.)
# For in put data use only proteinGroups.txt and evidence.txt





#=============discriminate_KR_peptides_v0-1-3.r======================

#! /Library/Frameworks/R.framework/Resources
#Enable command arguments
args=commandArgs(TRUE)
x=args[1]
message("Code written by CristianV.A.Munteanu C2019")

#####################################################################
# Originally developped at the Protein Analysis Facility 
# (University of Lausanne, Switzerland)
# Adapted at the Protein Chemistry Facility of IBRA - RO (ibar.ro)
# Please contact cristian.munteanu@biochim.ro for any supplementary 
# information.


# Input files
#
# evidence.txt
# lists all identified peptides in all identified forms 
# (with or without modifications)
#
# Version log
# Version 0.10: discriminate KR peptides, computes peptide median 
# for technical replicates, remove reverse and contaminants
# Version 0.11: changed script from col. indices to column names
# (can be used with evidence files with columns in other positions)
# Version 0.12: added new columns in the final output table, like id,
# peptide id, protein group id etc.
# Version 0.13: inserted command for export of deleted entries from
# the inputed file
# Version 0.14: added option for median translation of both raw and 
# normalized data
# Version 0.15: light channel data info added to the exported table
# and raw intensities
#  1: Sequence
#  2: Modifications
#  3: Leading.razor.protein
#  4. Gene.names
#  5. Protein.names
#  6: K.Count
#  7: R.Count
#  8. Leading.proteins
#  9: Peptide.ID
# 10: Protein.group.IDs
# 11: Ratio.H.M
# 12: Ratio.H.M.normalized
#
# Column description:
# 1: The identified AA sequence of the peptide.
#
# 2: Post-translational modifications contained within the identified
# peptide sequence.
#
# 3: The identifier of the best scoring protein, from the
# proteinGroups file this peptide is associated to.
#
# 4: Names of genes this peptide is associated with.
#
# 5: Names of proteins this peptide is associated with.
#
# 6: The number of instances of K contained within the sequence.
# The value for this can reliably be determined in the case of
# labeling partners based on the distance between the partners. These
# counts are used to solidify the peptide identification process.
#
# 7: The number of instances of R contained within the sequence.
# The value for this can reliably be determined in the case of
# labeling partners based on the distance between the partners.
# These counts are used to solidify the peptide identification
# process.
#
# 8: The identifiers of the proteins in the proteinGroups file, with
# this protein as best match, this particular peptide is associated
# with. When multiple matches are found here, the best scoring
# protein can be found in the 'Leading Razor Protein' column.
#
# 9: The identifier of the non-redundant peptide sequence.
#
# 10: The identifier of the protein-group this redundant peptide
# sequence is associated with, which can be used to look up the
# extended protein information in the file 'proteinGroups.txt'. As a
# single peptide can be linked to multiple proteins (e.g. in the
# case of razor-proteins), multiple id’s can be stored here separated
# by a semicolon. As a protein can be identified by multiple
# peptides, the same id can be found in different rows.
#
# 11: The ratio between two heavy and medium label partners.
#
# 12: Normalized ratio between two heavy and medium label partners. 
# The median of ratio sub-populations was shifted to 1.
#
# Exp. id: The id for the time points considered and their replicates.
# The time points are marked by the ‘T‘ letter and then by a number
# designating the exactly time point (e.g. 0,1,2,3 etc.). These are
# followed by letters under the format ‘_A‘,‘_B‘,‘_C‘ etc. which
# designates their technical replicates (injections).
#####################################################################

##########################
# Read input files
message("Reading file...")
x=read.delim(file =x,stringsAsFactors=F,)

message("Removing contaminants and reverse ids...")
write.table(x[which(x$Reverse=="+"|x$Potential.contaminant=="+"),],
"evEntriesDel.tab",sep="\t",row.names=F,col.names=T)
x=x[which(!x$Reverse=="+"&!x$Potential.contaminant=="+"),]

message("Filtering out peptides containing both K and R...")
write.table(x[which(!x$R.Count==0&!x$K.Count==0),],
"evEntriesDel.tab",sep="\t",row.names=F,col.names=F,append=T)
x=x[which(x$R.Count==0|x$K.Count==0),]

message("Filter out non-labeled peptides (terminal peptides)...")
write.table(x[which(x$R.Count==0&x$K.Count==0),],
"evEntriesDel.tab",sep="\t",row.names=F,col.names=F,append=T)
x=x[which(!x$R.Count==0|!x$K.Count==0),]

message("Log transforming ratios and removing redundant cols...")
x[,c("Ratio.H.M","Ratio.H.M.normalized","Ratio.M.L",
"Ratio.M.L.normalized","Ratio.H.L","Ratio.H.L.normalized",
"Intensity.L","Intensity.M","Intensity.H")]=
log(subset(x,select=c("Ratio.H.M","Ratio.H.M.normalized",
"Ratio.M.L","Ratio.M.L.normalized","Ratio.H.L",
"Ratio.H.L.normalized","Intensity.L","Intensity.M","Intensity.H")),2)
x=within(x,rm("Reverse","Potential.contaminant"))

# Split by labeled residue
Rtable=x[which(!x$R.Count==0&x$K.Count==0),]
Ktable=x[which(x$R.Count==0&!x$K.Count==0),]

# QC check of tables (should be TRUE)
# sum(nrow(Rtable),nrow(Ktable))==nrow(x)

# Create new column for each time point
Rtable$TimePoint=gsub("_.*","",Rtable$Experiment)
Ktable$TimePoint=gsub("_.*","",Ktable$Experiment)

message("Splitting by experiment label...")
Rlist=split(Rtable,Rtable$TimePoint)
Klist=split(Ktable,Ktable$TimePoint)

# Add labeled res. to list names
names(Rlist)=paste0("R",names(Rlist))
names(Klist)=paste0("K",names(Klist))

message("Aggregating technical repls. and extracting data for each 
time-point...")
combine=function(x){
 temp=aggregate(x[,c("Ratio.H.M","Ratio.H.M.normalized",
 "Ratio.M.L","Ratio.M.L.normalized","Ratio.H.L",
 "Ratio.H.L.normalized","Intensity.L","Intensity.M","Intensity.H")],
 x[,c("Sequence","Modifications","Leading.razor.protein",
 "Gene.names","Protein.names",
 "Leading.proteins","K.Count","R.Count","Protein.group.IDs",
 "Peptide.ID")],FUN=median,na.rm=T)
 return(temp)
}
Rlist=lapply(names(Rlist),function(x) combine(Rlist[[x]]))
names(Rlist)=gsub("K","R",names(Klist))
Klist=lapply(names(Klist),function(x) combine(Klist[[x]]))
names(Klist)=gsub("R","K",names(Rlist))

# Function to merge multiple dfs.
# Can also use merge with Reduce function from base R but needs
# to be adjusted
join=function(x){
 temp=merge(x[[1]],x[[2]],by=c("Sequence","Modifications",
 "Leading.razor.protein","Gene.names","Protein.names",
 "K.Count","R.Count","Leading.proteins","Peptide.ID",
 "Protein.group.IDs"),
 all=T,suffixes=paste0("_",names(x)[1:2]))
 for(i in 2:(length(x)-1)){
  temp=merge(temp,x[[i+1]],by=c("Sequence","Modifications",
  "Leading.razor.protein","Gene.names","Protein.names",
 "K.Count","R.Count","Leading.proteins","Peptide.ID",
 "Protein.group.IDs"),
  all=T,suffixes=paste0("_",names(x)[c(i,i+1)]))
  }
 return(temp)
}

Rdata=join(Rlist)
Kdata=join(Klist)

# Remove rows containing only NAs
Rdata=Rdata[!rowSums(is.na(Rdata[,11:64]))==54,]
Kdata=Kdata[!rowSums(is.na(Kdata[,11:64]))==54,]

# Center to median
Rdata[,which(grepl("Ratio",names(Rdata))&!grepl("ormalized",
names(Rdata)))]=scale(Rdata[,which(grepl("Ratio",names(Rdata))&!
grepl("ormalized",names(Rdata)))],scale=F,center=apply(Rdata[,which(
grepl("Ratio",names(Rdata))&!grepl("ormalized",names(Rdata)))],2,
median,na.rm=T))
Kdata[,which(grepl("Ratio",names(Kdata))&!grepl("ormalized",names(
Kdata)))]=scale(Kdata[,which(grepl("Ratio",names(Kdata))&!
grepl("ormalized",names(Kdata)))],scale=F,center=apply(Kdata[,which(
grepl("Ratio",names(Kdata))&!grepl("ormalized",names(Kdata)))],2,
median,na.rm=T))

Rdata[,which(grepl("Ratio",names(Rdata))&grepl("ormalized",
names(Rdata)))]=scale(Rdata[,which(grepl("Ratio",names(Rdata))&
grepl("ormalized",names(Rdata)))],scale=F,center=apply(Rdata[,which(
grepl("Ratio",names(Rdata))&grepl("ormalized",names(Rdata)))],2,
median,na.rm=T))
Kdata[,which(grepl("Ratio",names(Kdata))&grepl("ormalized",names(
Kdata)))]=scale(Kdata[,which(grepl("Ratio",names(Kdata))&
grepl("ormalized",names(Kdata)))],scale=F,center=apply(Kdata[,which(
grepl("Ratio",names(Kdata))&grepl("ormalized",names(Kdata)))],2,
median,na.rm=T))

message("Exporting tables...")
write.table(Rdata,
sep="\t",row.names=F,col.names=T,"R_pepsEvFile.tab")
write.table(Kdata,
sep="\t",row.names=F,col.names=T,"K_pepsEvFile.tab")

message("Finished")
