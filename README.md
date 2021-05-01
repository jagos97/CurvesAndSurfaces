



This assignment was finished on Visual studio, windows 10 with opengl 4.6. Bonus 1 was completed


This assignment has 2 different modes

Mode 1 is 2D:
	Here you can make bezier and b spline curves that can also be used for the 3D mode. You can click on the screen where there aren't any points to add a point.
	Clicking on a point you can select it and drag it around if you want. Once selected a point can also be deleted with 'X'.
	Using the right and left arrow keys you can change from bezier and b spline curves.
** Any number of control points can be added for bezier but only up to 43 points can be added for b-spline before it no longer works properly (will start messing up if you check out the parametric surface made with b spline).
** This because I couldn't quite make it work with variable sizes of the index array so I just used really big numbers. line 765 can be uncommented to display number
	of control points and indices used for making the curve and line 460 can be uncommented to check the number of indices used for the tensor and parametric surfaces 

Mode 2 is 3D
	This mode introduces a moving camera. By clicking and dragging you can rotate the camera to where you want to go. 'W' and "S" moves forwards and back relative to the direction you're looking.
	"A" and 'D' moves sideways relative to the position you're looking at and for added convience 'X' and 'C' can be used to move up and down as well.
	Pressing 'P' will change how you look at the curves (wirefram and fill).
	There are 4 scenes in this mode
		Scene 1 lets you see the curve and the control points in 3D
		Scene 2 lets yoou see the parametric surface created using the curve and control points created in 2D mode
		Scene 3 Gives you 3 tensor surfaces created using hard coded points using bezier curves. You can change which hard coded points to use using the up and down arrow keys
		Scene 4 gives you the same 3 tensor surfaces but using B-spline. 



Controls

W				-----------------------------		Move forward
S				-----------------------------		Move backwards
A				-----------------------------		Move Left
D				-----------------------------		Move Right

X				-----------------------------		Delete selected point	 /	Move up
C				-----------------------------		Move down

Space			-----------------------------		Change mode (2D, 3D)
P				-----------------------------		Change Polygon Mode 

Left-Click		-----------------------------		Add,Select or Move Point	/ Rotate camera
Left-Arrow		-----------------------------		Change Scene Decrement
Right-Arrow		-----------------------------		Change Scene Increment
Up-Arrow		-----------------------------		Change Tensor Surface Increment
Down-Arrow		-----------------------------		Change Tensor Surface Decrement
