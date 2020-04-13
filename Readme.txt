This is an implimentation of the consurf system
It takes a single chain pdb file called wild_mini.pdb (often provided from critires)
it does provides two ways to get consurf values one using the top 150 homologues 
and the other using the top 150 and then 150 other random matches
the default is to use the 300, 

To run you need a directory structure such as

/top/consurf_scripts
    /pdbidchain    e.g. 1fnfa
    /pdbidchain/wild_mini.pdb    - a single chain pdb file

you can then run it as
    ../consurf_scripts/consurf_home.sh

if you want both the 150 and 300 results you can then run
    ../consurf_scripts/consurf150.sh
if you only want the 150 results then edit consurf_home.sh to use those instead
