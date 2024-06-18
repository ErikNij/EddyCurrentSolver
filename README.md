# EddyCurrentSolver

This is the base case for induction heating simulation based off of [paper name]. While I did not write this code, I will attempt to make some sort of documenation to show how to use this code. 

# How to run

If you just want to run it, then you need to foam extend 4.1, and openfoam (.org or .com, any version (ish)) 

The tutorial "inductionHeating" under electromagnetics is what I have been working the most on. You can try to copy this to your run folder, and try following the steps below.

1)first source fe41, then do ./Allrun.coil, 
2)then you need to source openFoam to do the snappyHexMeshing. After sourcing, do ./Allrun.mesh, this can take a few minutes, but should not take more than 5. 
3)then you need to source fe41, for to do the rest of the simutions. The next few commands are done with ./Allrun <flag>. There are many flags, you can see them in the ./Allrun file if you want. They are:
clean
prepare
decompose
reset
reconstruct
start
force

4)so the first step is to clean the case, do this with ./Allrun clean
5)the next step is to prepare the case. This does some mesh manipulation, setting sets, and more, do this with ./Allrun prepare
6)the next step is decompose the case. do this with ./Allrun decompose
7)finally, we can now run the case, do this with ./Allrun start
8)to visulize the case, you can reconstruct, do this with ./Allrun reconstruct
i dont know what force and reset do at this moment. 

# Help installing 

This is the user folder of a foam-extend instalation. So first thing to do is set that up and complie everything. This probably wont work, beucase I have made some mistakes that I am working on fixing. 

There can be a lot of "bugs" / troubles while doing these steps, esspically for the first time, here are some tips / what I did . 

First, make sure that $FOAM_USER_LIBBIN is configured correctly, and all the $FOAM_USER_ cause when it does the make lists, it might not be correct
In steps 3) if you get the error something like expected faceList got comapactFaceList, then you need to make a function in normal OpenFOAM to convert between the 2, you can find it online.
On step 5), if the first thing fails imdeitly then you didnt compile the fe41 user part correctly, see the first tip, and make sure you made it in the first place
If step 6) fails with the error cannot find something.so, then you need to make sure the libbin is correct, and add the .so file in to your control dict

# Editing the base case

If you want to change things, here is where you can do that.
If you want to  change the size of the (background) mesh this is done @ constant/polyMesh.org/blockMeshDict.py. Make sure you are eding the .org folder
To change the coil, go to constant/featureEdgeMesh.org/coil.py. Here you can make the coil wider or have more turns. If you add more turns then you also need to change the file at constant/edgeBiotSavartProperties
To change the frequency or current ampliture, go to f0 and edgeBiotSavartProperties respectivly. 

To change the work piece, follow the following instrictions. The best way to make a mesh is in blender, make the wp, and then get the inverse. The inverse is not strictly nesscaarry, but it is used to make a more refined region around the WP.

You may be wondering how to create an stl file, and esspeically the inverse. Well, here is how. 

Lets say you have a workpeice, thats good. lets assume its in stl. Go to blender, open that stl. 
then add a cube, this can be done by going to modeling > add > mesh > cube
you will likly need to resize it. Currently, sizes are x = 104mm, y = 154mm z = 135mm. This can be done by going to the right to the orange box, and change the scale. However, the scale is half the dimetions of the cube, so you need:
X:0.052
Y:0.077
Z:0.0675

You will also need to move it in the postive z direction y 0.0225.

to do the inverse, you will need the boolean modifer. To add a modifer you go to the modifers tab on blender, which the wrench one below the orange box.
the add modifer > boolean. Then make sure difference is selected. For object, select your workpeice. I dont know what solver is better, but I did exact. It should then be done.
You can then export the buffer as an stl. You need to make sure that the workpeice is not exported with it. Thselect the buffer, and then while exporting, check "export only selection".
Done!

# Tips, tricks, and runtime debugging

If the solver fails, here are some things to check.

1) If it fails around the trying to calculate the A0, then your coil could be illdefined. Make sure that the coil r is slightly greater than x & y. 
