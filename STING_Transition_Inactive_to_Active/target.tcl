mol addfile active.pdb 
set all [atomselect top "all"]
$all set beta 0
set prot [atomselect top "all "]
$prot set beta 1
$all writepdb target.cnst
exit

