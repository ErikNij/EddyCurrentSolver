# EddyCurrentSolver

This is the base case for induction heating simulation based off of [paper name]. While I did not write this code, I will attempt to make some sort of documentation to show how to use this code. 

# How to run

If you just want to run it, then you need to foam extend 4.1, and OpenFOAM (.org or .com, any version (ish)) 

The tutorial "inductionHeating" under electromagnetics is what I have been working on the most on. You can try to copy this to your run folder and follow the steps below.

1. To run the most simple case, we need to prepare the coil, mesh, and sets, then run the case. To do this, copy the tutorial to your run folder, and navigate to the folder. First source `fe41`, then do `./Allrun.coil`, now the coil is created.
2. Now we need to mesh the case. To do this, we use snappyHexMesh, however foam-extend does not have the most up to date version of snappyHexMesh, so we will switch to OpenFOAM (.org or .com) to use it. to do this, you need to source `of` After sourcing, do `./Allrun.mesh`. This can take a few minutes, depending on the complexity of the mesh you have and the detail you are trying to render it at.  
3. For the rest of the steps, we will use foam-extend and other custom functions that are saved in the foam-extend folder, so source `fe41` for to do the rest of the steps. The next few commands are done with `./Allrun <flag>`. There are many flags, you can see them in the ./Allrun file. They are:
   - `clean`
   - `prepare`
   - `decompose`
   - `reset`
   - `reconstruct`
   - `start`
   - `force`
5. The first step is to clean the case. Do this with `./Allrun clean`
6. The next step is to prepare the case. This does some mesh manipulation, setting sets, and more, do this with `./Allrun prepare`
7. The next step is decomposing the case. Do this with `./Allrun decompose`
   - Note: This step is not currently used
9. Finally, we can now run the case. Do this with `./Allrun start`
10. To visualize the case, you can reconstruct. Do this with `./Allrun reconstruct`
    - Note: This step is also currently not used. 
I don't know what force and reset do at this moment.

Furthermore, all of these commands are combined in one bash file called setup.sh, which can be run with `bash setup.sh`

# Help installing 

This is the user folder of a foam-extend installation. So the first thing to do is set that up and compile everything. This probably won't work, because I have made some mistakes that I am working on fixing. 

There can be a lot of "bugs" / troubles while doing these steps, especially for the first time. Here are some tips / what I did. 

First, make sure that $FOAM_USER_LIBBIN is configured correctly, and all the $FOAM_USER_ cause when it does the make lists, it might not be correct
In steps 3) if you get the error something like expected faceList got comapactFaceList, then you need to make a function in normal OpenFOAM to convert between the 2, you can find it online.
On step 5), if the first thing fails imdeitly then you didn't compile the fe41 user part correctly, see the first tip, and make sure you made it in the first place
If step 6) fails with the error cannot find something.so, then you need to make sure the libbin is correct, and add the .so file in to your control dict

# Editing the base case

If you want to change things, here is where you can do that.
If you want to  change the size of the (background) mesh, this is done at constant/polyMesh.org/blockMeshDict.py. Make sure you are editing the .org folder
To change the coil, go to constant/featureEdgeMesh.org/coil.py. Here, you can make the coil wider or have more turns. If you add more turns, then you also need to change the file at constant/edgeBiotSavartProperties
To change the frequency or current amplitude, go to f0 and edgeBiotSavartProperties respectively. 

To change the work piece, follow the following instructions. The best way to make a mesh is in blender, make the wp, and then get the inverse. The inverse is not strictly necessary, but it is used to make a more refined region around the WP.

You may be wondering how to create an STL file,  especially, the inverse. Well, here is how. 

Let's say you have a workpiece, that's good. Let's assume it's in STL. Go to blender and open that STL. 
then add a cube, this can be done by going to modeling > add > mesh > cube
you will likely need to resize it. Currently, sizes are x = 104mm, y = 154mm z = 135mm. This can be done by going to the right to the orange box, and change the scale. However, the scale is half the dimensions of the cube, so you need:
`X:0.052`
`Y:0.077`
`Z:0.0675`

You will also need to move it in the positive z direction by `0.0225`

To do the inverse, you will need the boolean modifier. To add a modifier you go to the modifiers tab on blender, which is the wrench one below the orange box.
The add modifier > boolean. Then make sure the difference is selected. For object, select your workpiece. I don't know what solver is better, but I did exact. It should then be done.
You can then export the buffer as an STL. You need to make sure that the workpeice is not exported with it. Then select the buffer, and then while exporting, check "export only selection".
Done!

# Tips, tricks, and runtime debugging

If the solver fails, here are some things to check.

1) If it fails around the trying to calculate the A0, then your coil could be ill defined. Make sure that the coil r is slightly greater than x & y.
2) If the snappyHexMesh fails saying:
   `--> FOAM FATAL ERROR: (openfoam-2212)
   Shell geometry_buffer.stl does not support testing for inside
   Probably it is not closed.`

   Then, you either have an open mesh, or overlapping meshs. Open meshs can be fixed by looking up online, however this is rarly the problem. Often, your mesh have overlapping sections. To fix this, load up the stl of your workpeice, and in blender choose the remesh modifer under generate. Keep the mode on voxel, and decrease the voxel size until you are happy with the result. Note: make sure to do this before doing the boolean with the buffer, as remeshing while having the boolean on takes a lot more computing power.
3) If your mesh is not smooth, then it could be a problem with the set mesh region. Instead of only looking at the conductor, look at the entire internal mesh, and see if there is a smooth line that you want somewhere. If you do see it somewhere, then it is a problem with setSetBatch. This sets the regions that are used later. If you open setSetBatch in the system folder, you will see the pointSet command, this is where points get assigned to a region. All you have to do it set the tolerance to be lower. The tolerance should be a bit lowerer then the size of a cell at the edge (in the r direction). Note: the tolerance can also be set too low, it has to just be a bit lower than the cell size.
