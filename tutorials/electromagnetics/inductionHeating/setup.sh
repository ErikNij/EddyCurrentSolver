#! /bin/sh.
source /home/erik/foam/foam-extend-4.1/etc/bashrc
./Allclean
./Allrun.coil
source /opt/OpenFOAM/OpenFOAM-v2212/etc/bashrc
./Allrun.mesh
source /home/erik/foam/foam-extend-4.1/etc/bashrc
./Allrun prepare
