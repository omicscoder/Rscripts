MS/MS spectra analysis of interpeptide disulfide bridge analysis; 

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
