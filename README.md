# Script gen_insilico_disulfide_peps_v03.r

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



# Script calcProtFDR_v0-1-2.r

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


# Script calcpI.r

This script reads amino acid sequences of peptides or proteins and
returns their theoretical pI

The returned values does NOT take into account co-translational or
post-translational modifications (neither chemical modifications)


# Script discriminate_KR_peptide_v0-1-3.r

pcSILAC data analysis of A375-ST cells labeled with R6/R10,K0/K0 and chased for various 
times (t0-t5) in R0/R0,K4/K8 after EDEM2 induction with TET (see DOI: 10.1016/j.mcpro.2021.100125)

The data analysis workflow will closelly follow the one described 
Munteanu CVA. et. al., Munteanu CVA#, Chiritoiu GN#, Chiritoiu M, Ghenea S, Petrescu AJ, Petrescu SM., 
"Affinity Proteomics and Deglycoproteomics Uncover Novel EDEM2 Endogenous Substrates and an Integrative ERAD Network.",
Mol. Cell. Proteomics; 20: 100125, 2021

For in put data use only proteinGroups.txt and evidence.txt



