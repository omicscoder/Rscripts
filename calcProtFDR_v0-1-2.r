# Code for the script calcProtFDR_v0-1-2.r
#
# This script reads a txt (tab separated data) with a list of protein
# and exports a short list at 1% or 5% FDR estimated at the protein 
# level
# The protein list should be obtained after searching both, a target
# and a decoy database
# The first column of the table should contain the 
# UniProtID/Accession number
# Decoy entries should have the 'DECOY'tag in the first  column 
# Protein FDR is calculated as described in Reidegeld et al., 
# Proteomics. 2008 Mar;8(6):1129-37




#=============calcProtFDR_v0-1-2.r======================

#! /Library/Frameworks/R.framework/Resources
#Enable command arguments
args=commandArgs(TRUE)
x=args[1]
y=args[2]
message('Code written by CristianV.A.Munteanu C2019')
message('For FDR methodology see Reidegeld et al.,Proteomics. 2008 
Mar;8(6):1129-37')
######################################################################
# Developped at the Protein Protein Chemistry Facility of IBRA - RO
# (https://www.biochim.ro/facility-1/) by Dr. Cristian (V.A.) Munteanu
# Please contact cristian.munteanu@biochim.ro
# (http://cristianvamunteanu.tk) for any supplementary information.


# Input files
#
# ProteinGroupList.txt/FDR_correctedProteinGroupList.tab
# results to be exported
#
# Version log
# Version 0.10: reads in a tab formatted list with protein groups 
# usually exported from Proteome Discoverer and returns a shortlist of
# of protein groups by applying a protein FDR (pFDR) as described in 
# the Reidegeld et al. journal article
# Version 0.11: added option to analyze files containing individual 
# scores from multiple searches
# Version 0.12: export option was modified to keep the original input
# filename to which the suffix '_scoreFDR' is appended
######################################################################

##########################
t1=Sys.time()
message("Loading functions...")

# Define the function for pFDR calculation
pFDR=function(x){
# Create matrix to store res.
 temp=matrix(,length(x),2)
# Populate the table with pFDR scores
 for(i in length(x):1){
  temp[i,1]=i
  temp[i,2]=length(x[which(grepl('DECOY',x[1:i]))])/length(x[1:i])
  }
# Return the results
 return(temp)
}

##########################
# Read input files
message('Reading input file...')
res=read.delim(file =x,stringsAsFactors=F,)

# Calculate current protein FDR
message('Calculating protein FDR...')
message(paste0('Current protein FDR is ',
round(length(res[which(grepl('DECOY',
res[,1])),1])/length(res[,1])*100,2),'%'))

# Calculate protein identification maximum score
message('Calculating protein scores...')
res$MaxScore=if(length(res[which(grepl('core',names(res)))])==1)
{res[,which(grepl('core',names(res)))]}else{
apply(X=res[,which(grepl('core',names(res)))],1,FUN=function(x)
max(x,na.rm=T))}

# Sort protein groups by score
res=res[order(res$MaxScore,decreasing=T),]

# Get pFDR distribution
temp=pFDR(res[,1])

# Subset data according to the selected pFDR
res=res[1:tail(temp[which(temp[,2]<=y),],1)[,1],]

# Filter only for target db. matches
res=res[which(!grepl('DECOY',res[,1])),]

# Export the results
message("Writing results tables...")
write.table(res[,-ncol(res)],
paste0(gsub('.txt','',x),'_',gsub('\\.','',y),'pFDR','.tab')
,sep='\t',row.names=F,col.names=T)
message(paste0('New protein FDR is ',
round(tail(temp[which(temp[,2]<=y),],1)[,2]*100,2),'%'))

message("Finished!")
t2=Sys.time()
message(paste0('Protein list processed in ',round(t2-t1,2),' secs.'))
##########################
