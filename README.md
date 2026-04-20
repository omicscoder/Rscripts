gen_insilico_disulfide_peps_v03.r

The script can take as input a list of tryptic peptides and obtain in silico
the MS/MS tables with all theoretical b adn y ions of pairs between two
tryptic peptides.

This script reads a txt/tab file containing a list of tryptic peptides and 
exports the the PeptideTheoreticalIonsTable from the two peptide pairs;
the script is based on Xu H, Zhang L, Freitas MA., J Proteome Res. 2008,
Jan;7(1):138-44. doi: 10.1021/pr070363z.

Requirements: 'OrgMassSpecR' package

Plase install the package using 'install.package()' function if this is 
not yet available in your R installation

Script calcProtFDR_v0-1-2.r
This script reads a txt (tab separated data) with a list of protein
and exports a short list at 1% or 5% FDR estimated at the protein 
level
The protein list should be obtained after searching both, a target
and a decoy database
The first column of the table should contain the 
UniProtID/Accession number
Decoy entries should have the 'DECOY'tag in the first  column 
Protein FDR is calculated as described in Reidegeld et al., 
Proteomics. 2008 Mar;8(6):1129-37

