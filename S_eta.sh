#!/bin/bash
#spherical density -> S_eta
#this code should run after finishing spherical density analysis 

lvalue=6  #mode of spherical harmonics to be analyzed; 3 for silica and 6 for bljm
corr=13  #correlation atyp1_atyp2
frac=0.01  #fraction of the top and bottom cells to define rho_max and rho_min; a treatment to maintain the stability of the two extrema 

path=./structure/   #path to load the spherical density files
cd $path
rm -f S_eta-vs-r-l$lvalue-$corr-$frac.dat  #delete the old file if exist

#loop over a serie of density files for the S_eta analysis
for file in $(ls spherical-density-r*-bljm-corr$corr.xyz); do  
 rad=${file:19:5}  #r read from the file name
 sleep 1s

 awk 'NR>2{print $2,$3,$4,$5}' $file >tmpfile.dat  #x,y,z,rho
 nlines=`less tmpfile.dat |wc -l`

 #rho_max and rho_min as the mean value of the top and bottom 1% cells
 maxval=`sort -grk 4 tmpfile.dat |awk 'NR<='$nlines'*'$frac' { sum1 += $4; n1++ } END { if (n1 > 0) printf("%f\n", sum1 / n1); }'`
 minval=`sort -gk 4 tmpfile.dat |awk 'NR<='$nlines'*'$frac' { sum2 += $4; n2++ } END { if (n2 > 0) printf("%f\n", sum2 / n2); }'`

 #calculate the normalized density field (only count cells with positive density values)
 awk '{($4=($4-'$minval')/('$maxval'-'$minval'));if($4>0) print ;}' tmpfile.dat >normalized-density.dat

 #cp tmpfile.dat normalized-density.dat  #the orginal spherical density, should yield S_rho 
 ncells=`less normalized-density.dat |wc -l`  
 echo max_rho $maxval, min_rho $minval, n_cells $ncells

 #run the fortran code for obtaining S_eta
 gfortran -o seta.out ../S_eta.f90 && ./seta.out $lvalue $ncells
 s_eta=`awk '{print}' S_eta.dat`

 #output S_eta vs. r 
 echo $rad $s_eta $file
 echo $rad $s_eta >> S_eta-vs-r-l$lvalue-$corr-$frac.dat
done


