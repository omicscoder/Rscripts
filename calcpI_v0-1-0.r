# Code for the script calcpI.r
# 
# This script reads amino acid sequences of peptides or proteins and
# returns their theoretical pI
# The returned values does NOT take into account co-translational or
# post-translational modifications (neither chemical modifications)
#
#



#=======================calcpIR_v0-1-0.r==============================

#! /Library/Frameworks/R.framework/Resources
#Enable command arguments
args=commandArgs(TRUE)
x=args[1]
#y=args[2]
message('Code written by CristianV.A.Munteanu C2019')
######################################################################
# Developped at the Protein Protein Chemistry Facility of IBRA - RO
# (https://www.biochim.ro/facility-1/) by Dr. Cristian (V.A.) Munteanu
# Please contact cristian.munteanu@biochim.ro
# (http://cristianvamunteanu.tk) for any supplementary information.


# Input files
#
# FASTAseqs.txt/seqs.tab
#
# Version log
# Version 0.10: reads in a FASTA formatted list with protein seqs. 
# and returns a table with calculated pI for each entry (default to
# Bjellqvist pK values)
# Version 0.11: should allow selecting the pK scale similar with the
# ones available in 'Peptides' package
######################################################################

##########################
t1=Sys.time()
message("Loading functions...")

# Load aa. pK data (Bjellqvist scale)
aapK=data.frame(cbind(
aa=as.character(c("C","D","E","H","K","R","Y","cTer","nTer")),
pK=as.numeric(c(9.00,4.05,4.45,5.98,10.00,12.00,10.00,3.55,7.50))
),stringsAsFactors=F)
aapK$pK=as.numeric(aapK$pK)

# Function to count no. of charged aa. occurences
countaa=function(x){
# 'x' - peptide/protein aaseq.
 temp=as.numeric()
 for(i in 1:7){
  temp[i]=nchar(x)-nchar(gsub(aapK$aa[i],'',x))
 }
 return(temp)
}
countaa('ALVCC')

# Function to calculate pI of hydrophobic aa. containing peps.
# 'x' - peptide/protein aas seq
pI=function(x){
# Create numeric vector to store the count of each  charged aa.
temp=as.numeric()
# Count the frequency of each aa.
for(i in 1:7){
  temp[i]=nchar(x)-nchar(gsub(aapK$aa[i],'',toupper(x)))
}
# Calculate pI
as.numeric(uniroot(
function(x) 
# N-terminus
1/(1+10^(x-7.5))+
# Cys
temp[1]*(-1/(1+10^(9.00-x)))+
# Asp
temp[2]*(-1/(1+10^(4.05-x)))+
# Glu
temp[3]*(-1/(1+10^(4.45-x)))+
# His
temp[4]*1/(1+10^(x-5.98))+
# Lys
temp[5]*1/(1+10^(x-10.00))+
# Arg
temp[6]*1/(1+10^(x-12.00))+
# Tyr
temp[7]*(-1/(1+10^(10.00-x)))+
# C-terminus
-1/(1+10^(3.55-x))
,c(0,14))[1])
}


# Define the function for pI calculation
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
