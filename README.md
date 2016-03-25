# Graphics
Projects done for COEN 290 - Computer Graphics

Project#1:
Implement an OpenGL graphics program with the following capabilities: 1. Creates a computer-generated 3D character object (e.g., robots or insects with head,
body, arms, and legs). 2. Renders the 3D character in a 3D scene. 3. Displays the 3D character in wireframe representation. 4. Supports the following keyboard inputs:
a. ‘a’ – apply the subsequent parameter updates to the ‘At’ point
b. ‘f’ – apply the subsequent parameter updates to the ‘From’ point  
c. ‘x’ – apply the subsequent increments/decrements to the x coordinate value 
d. ‘y’ – apply the subsequent increments/decrements to the y coordinate value 
e. ‘z’ – apply the subsequent increments/decrements to the z coordinate value 
f. ‘i’ – increase the current x, y or z coordinate by delta 
g. ‘d’ – decrease the current x, y or z coordinate by delta 
h. ‘r’ – reset the At point to the default position 
i. ‘c’ – reset the From point to the default position

Project#2:
Implement an OpenGL graphics program with the following capabilities:
1. Use hierarchical modeling to ‘connect’ all the parts of your 3D character object. 
2. Allow the user to select any part of your 3D character object (via pop-up menu, keyboard input, or selection/picking) for subsequent transformations. 
3. Use OpenGL lighting to shade the 3D character object. 
4. Allow the user to record/playback user input events by accepting the following entries from a pop-up menu. 
  a. Start/begin recording: store the subsequent transformation events 
  b. Stop/end recording 
  c. Load/read playback_<your initials>.txt, which contains user input events from a previous session. 
  d. Playback the saved events 
  e. Reset the default positions/orientations of your 3D character object in the scene 
  f. Save the current recording to playback_<your initials>.txt
5. Create an interesting animation (at least 20 seconds).

Project#3:
Write an interactive curve design program that allows the user to draw B-spline curves. Your program will provide a pop-up menu with the following commands:
1. Control point input On/Off.  When the control point input is ON, the user will be allowed
to add control points. 
2. Control polygon On/Off. When this option is ON, your program will display the current control polygon. 
3. B-spline curve On/Off. When this option is ON, your program will display the B-spline curve based on the current control points. 
4. Select a control point: In this mode, your program will allow the user to select a control point.   
5. Delete or move the selected control point.  Your program will automatically update the current control polygon/B-spline curve 
when the selected control point is deleted or moved. 
6. Save: Save the current control points to bspline.txt 
7. Retrieve: Retrieve control points from bspline.txt 
8. Clear: Clear the current display window and delete all control points. 
9. Draw wireframe surface: Display the B-spline surface as a wireframe mesh. 
10. Shade surface: Shade the B-spline surface using pre-specified lighting and material parameters. 
11. Texture surface: Texture map an image onto the B-spline surface.

Project#4
Implement a recursive ray tracing program to generate an interesting scene.  In addition to spheres and planes, your program needs
to include at least one surface mesh generated from Project #3. Implement the Phong lighting model for computing the basic object 
color. Your ray tracer will support reflection and shadow.   The minimum depth of the recursion is 3.
