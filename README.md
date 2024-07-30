# 3D Mesh Subdivision and Visualization

## Overview
This project is an exploration of surface meshes and their modifications through subdivision. Using Processing and Java, the program reads polyhedral models and creates Loop subdivided surfaces. It also demonstrates the "corners" representation by visualizing corner operations such as next, opposite, and swing.

### Features
- Mesh Loading: Load different mesh files (tetrahedron, octahedron, icosahedron, star, torus) with a keystroke.
- Loop Subdivision: Perform Loop subdivision on the mesh to create smoother surfaces.
- Corner Visualization: Visualize corner operations using small spheres to indicate the current corner.
- Shading Modes: Switch between flat shading and Gouraud shading using per-vertex normals.
- Interactive Rotation: Rotate the mesh by clicking and dragging the mouse.
### Loop Subdivision
Loop subdivision is a method for refining triangle meshes to produce smoother surfaces. Each time you subdivide the mesh, new vertices are created for each edge, and the positions of both new and old vertices are recalculated. The original triangles are replaced by smaller triangles, enhancing the mesh's smoothness.

### Subdivision Steps
1. Create New Vertices: New vertices are added at the midpoint of each edge, with their positions calculated based on surrounding vertices.
2. Update Old Vertices: The positions of the original vertices are adjusted to maintain the mesh's overall shape and smoothness.
3. Reconstruct Triangles: Each original triangle is divided into smaller triangles using the new vertices.
## Key Commands
- 1-5: Load different mesh files (tetrahedron, octahedron, icosahedron, star, torus).
- d: Perform Loop subdivision on the current mesh (can be done multiple times).
- f: Display the mesh with flat shading.
- g: Display the mesh with Gouraud shading using per-vertex normals.
- n: Change the current corner using the "next" operator.
- p: Change the current corner using the "previous" operator.
- o: Change the current corner using the "opposite" operator.
- s: Change the current corner using the "swing" operator.
- q: Quit the program.
### How to Run
1. Download and install Processing.
2. Open the provided .pde file in the Processing IDE.
3. Run the sketch to load and manipulate 3D meshes.
