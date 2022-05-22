# three-dimensional-structure

this fortran codes perform the 3d structural order analysis for a binary LJ mixture (bljm) and 
silica-based network system (nsx) as documented in the paper of Zhang and Kob (see ref at the end). 
The configuration (in xyz format) to be analyzed is by default 3d-periodic, but non-periodic 
configs can easily be adapted to this analysis. 

particle types: 1->big; 2->small (bljm)  1->O; 2->Si; 3->Na (nsx)

usage of the code:

gfortran -o bljm.out 3d-structure-analysis-bljm.f90 && ./bljm.out sff1 11 1.5 2.5 0.05 0.2 2000 6 0.1

As one can see from above, runing the code requires to load some parameters, which define which/how the 
analysis should be performed. see the related segement in the code for their definition and recommended values

local coordinate system: 

this lies at the core of the proposed four-point correlation function. the general rule 
is that one takes a relatively stable local structure motif, formed by a triplet of atoms,
as reference. as such, one first needs some basical knowledge of the structure of 
the system under investigation, such as g(r).
specifically, for hard-shere-like systems, one can use the triplet of particles that are 
nearest neighbor to each other to construct the reference coordinate system; 
for open-network systems with a well-defined local structure, e.g. SiO4 in sio2-based glasses,
one can use the O-Si-O linkage to construct the local reference frame
in any case, one should play a bit around with the choices of local reference frames
to see their performance

output of the analysis:

1) spherical density distribution (aka sff2 in the code): 
project the spatial distribution of the particles on a sphere with respect to
the chosen local coordinate system. the code will write a xyz file for visualization
of the 3d density distribution, using tools such as OVITO.

2) quantitative characterization of the 3d structure (aka sff1 in the code):
decompose the signal as seen from the density distribution by spherical harmonics (SH), 
for which different modes (l) catch different types of symmetries and order. the code 
will write a dat file of r vs. S_rho.

Reference: Zhang, Z., & Kob, W. (2020). Revealing the three-dimensional structure of liquids
using four-point correlation functions. PNAS, 117(25), 14032-14037.

for more information concerning the usage of the code and the proposed four-point correlation
function, feel free to contact Zhen Zhang by email: zhen.zhang1991@hotmail.com
