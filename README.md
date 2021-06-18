# Hard-fit-by-minimization
An efficient Fortran program to fit a stress-strain curve by the empirical function derived from crystal plasticity.

This is a Fortran program to fit an experimental strain hardening curve by a phenomenological  program proposed in a paper submitted to Materials Sceince and Enginnering A:

Title of the paper:
A new macroscopic strain hardening function based on microscale crystal plasticity and its application in polycrystal modeling

Authors and affiliation:
Sudeep K. Sahoo,a,b, Satyaveer Singh Dhinwal,a,b, Viet Q. Vu,a,b,c, Laszlo S. Toth,a,b,d

a Laboratoire d’Etude des Microstructures et de Mécanique des Matériaux (LEM3), CNRS UMR 7239, Université de Lorraine, 57045 Metz Cedex 1, France

b Laboratory of Excellence on Design of Alloy Metals for low-mAss Structures (DAMAS), Université de Lorraine, Metz, France

c Thai Nguyen University of Technology, Thai Nguyen, Vietnam

d University of Miskolc, Institute of Physical Metallurgy, Metal Forming and Nanotechnology, Miskolc, Hungary.


First provide an ASCI file containing 10-50 experimental points taken from your measured strain hardening curve. 

An example is given in the file: exp-curve.dat.

Then edit the hard.ctl file for your expected parameter range and options. 

Run the program hard.exe.

There will be thre output files:

Result.dat: Results obtained as a function of strain: stresses and derivative for Kocks-Mecking plot

Hardplot.dat: The simulated stress-strain curve and the Kocks-Mecking curve ready for plotting

Running.dat: Information on the running, it contains also the values of teh parameters
