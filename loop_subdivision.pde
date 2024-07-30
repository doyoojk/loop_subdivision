// Author: Jamie Kim
import java.util.Arrays;

// parameters used for object rotation by mouse
float mouseX_old = 0;
float mouseY_old = 0;
PMatrix3D rot_mat = new PMatrix3D();
boolean start, perVertex, gShading;
Mesh mesh;

// initialize stuff
void setup() {
  size (800, 800, OPENGL);
}

// Draw the scene
void draw() {
  
  background (170, 170, 255);
  
  // set up for perspective projection
  perspective (PI * 0.333, 1.0, 0.01, 1000.0);
  
  // place the camera in the scene
  camera (0.0, 0.0, 5.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);
    
  // create an ambient light source
  ambientLight (100, 100, 100);
  
  // create two directional light sources
  lightSpecular (204, 204, 204);
  directionalLight (102, 102, 102, -0.7, -0.7, -1);
  directionalLight (152, 152, 152, 0, 0, -1);
  
  pushMatrix();
  
  applyMatrix (rot_mat);  // rotate object using the global rotation matrix
  
  ambient (200, 200, 200);
  specular(0, 0, 0);
  shininess(1.0);
  
  stroke (0, 0, 0);
  fill(200, 200, 200);

  // THIS IS WHERE YOU SHOULD DRAW THE MESH
  
  // Display the square
  if (!start) {
    beginShape();
    normal(0.0, 0.0, 1.0);
    vertex(-1.0, -1.0, 0.0);
    vertex( 1.0, -1.0, 0.0);
    vertex( 1.0,  1.0, 0.0);
    vertex(-1.0,  1.0, 0.0);
    endShape(CLOSE);
  } else {
    drawMesh(mesh); 
    visualizeCorner(mesh.currCorner);
  }
    
  popMatrix();
}

void drawMesh(Mesh mesh) {
  for (int i = 0; i < mesh.V.size(); i += 3) {
    beginShape(TRIANGLES);
    for (int j = 0; j < 3; j++) {
      if (gShading) {
        noStroke();
        Vertex vertexNormal = mesh.vertexNorms[mesh.V.get(i + j)];
        normal(vertexNormal.x, vertexNormal.y, vertexNormal.z);
      } else {
        // Flat shading 
        Vertex normalFace = mesh.faceNorms.get(i / 3);
        normal(normalFace.x, normalFace.y, normalFace.z);
      }
      Vertex vert = mesh.G.get(mesh.V.get(i + j));
      vertex(vert.x, vert.y, vert.z);
    }
    endShape(CLOSE);
  }
}

void visualizeCorner(int corner) {
  int faceStart = (corner / 3) * 3;
  int vertexIndex1 = mesh.V.get(faceStart);
  int vertexIndex2 = mesh.V.get(faceStart + 1);
  int vertexIndex3 = mesh.V.get(faceStart + 2);

  // vertices of face
  Vertex v1 = mesh.G.get(vertexIndex1);
  Vertex v2 = mesh.G.get(vertexIndex2);
  Vertex v3 = mesh.G.get(vertexIndex3);

  // Calculate centroid of face with currCorner
  float centroidX = (v1.x + v2.x + v3.x) / 3;
  float centroidY = (v1.y + v2.y + v3.y) / 3;
  float centroidZ = (v1.z + v2.z + v3.z) / 3;
  
  // Get the position of the current corner vertex
  Vertex cornerVertex = mesh.G.get(mesh.V.get(corner));

  // Calculate the offset position for the sphere
  float sphereX = centroidX + (cornerVertex.x - centroidX) * 0.8;
  float sphereY = centroidY + (cornerVertex.y - centroidY) * 0.8;
  float sphereZ = centroidZ + (cornerVertex.z - centroidZ) * 0.8;

  // Draw the sphere
  pushMatrix();
  translate(sphereX, sphereY, sphereZ);
  fill(255, 0, 0);
  noStroke();
  sphere(0.04);
  popMatrix();
}

// handle keyboard input
void keyPressed() {
  if (key == '1') {
    read_mesh ("tetra.ply", 1.5);
    start = true;
  }
  else if (key == '2') {
    read_mesh ("octa.ply", 2.5);
    start = true;
  }
  else if (key == '3') {
    read_mesh ("icos.ply", 2.5);
    start = true;
  }
  else if (key == '4') {
    read_mesh ("star.ply", 1.0);
    start = true;
  }
  else if (key == '5') {
    read_mesh ("torus.ply", 1.6);
    start = true;
  }
  else if (key == 'n') {
    // next corner operation
    mesh.currCorner = mesh.NextCorner(mesh.currCorner);
  }
  else if (key == 'p') {
    // previous corner operation
    mesh.currCorner = mesh.PrevCorner(mesh.currCorner);
  }
  else if (key == 'o') {
    // opposite corner operation
    mesh.currCorner = mesh.OppositeCorner(mesh.currCorner);
  }
  else if (key == 's') {
    // swing corner operation
    mesh.currCorner = mesh.SwingCorner(mesh.currCorner);
  }
  else if (key == 'f') {
    // flat shading, with black edges
    gShading = false;
  }
  else if (key == 'g') {
    // Gouraud shading with per-vertex normals
    gShading = true;
  }
  else if (key == 'd') {
    // subdivide mesh
    mesh = mesh.subdivide(mesh);
  }
  else if (key == 'q') {
    // quit program
    exit();
  }
}

// Read polygon mesh from .ply file
//
// You should modify this routine to store all of the mesh data
// into a mesh data structure instead of printing it to the screen.
void read_mesh (String filename, float scale_value)
{
  String[] words;

  String lines[] = loadStrings(filename);

  words = split (lines[0], " ");
  int num_vertices = int(words[1]);
  println ("number of vertices = " + num_vertices);

  words = split (lines[1], " ");
  int num_faces = int(words[1]);
  println ("number of faces = " + num_faces);
  
  mesh = new Mesh(num_vertices, num_faces);

  // read in the vertices
  for (int i = 0; i < num_vertices; i++) {
    words = split (lines[i+2], " ");
    float x = float(words[0]) * scale_value;
    float y = float(words[1]) * scale_value;
    float z = float(words[2]) * scale_value;
    println ("vertex = " + x + " " + y + " " + z);
    mesh.G.add(new Vertex(x,y,z));
  }

  // read in the faces
  for (int i = 0; i < num_faces; i++) {
    
    int j = i + num_vertices + 2;
    words = split (lines[j], " ");
    
    int nverts = int(words[0]);
    if (nverts != 3) {
      println ("error: this face is not a triangle.");
      exit();
    }
    
    int index1 = int(words[1]);
    int index2 = int(words[2]);
    int index3 = int(words[3]);
    println ("face = " + index1 + " " + index2 + " " + index3);
    
    mesh.V.add(index1);
    mesh.V.add(index2);
    mesh.V.add(index3);
  }
  // Construct data tables
  mesh.OTable();
  mesh.CalcFaceNorms();
  mesh.CalcVertexNorms();
}

class Vertex {
  float x, y, z;

  Vertex(float x, float y, float z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  Vertex sub(Vertex second) {
    return new Vertex(x - second.x, y - second.y, z - second.z);
  }

  Vertex add(Vertex second) {
    return new Vertex(x + second.x, y + second.y, z + second.z);
  }

  Vertex cross(Vertex second) {
    return new Vertex(y * second.z - z * second.y, z * second.x - x * second.z, x * second.y - y * second.x);
  }

  Vertex norm() {
    float length = sqrt(x * x + y * y + z * z);
    return (length != 0) ? new Vertex(x / length, y / length, z / length) : new Vertex(0, 0, 0);
  }

  Vertex mult(float scale) {
    return new Vertex(x * scale, y * scale, z * scale);
  }

  boolean equals(Vertex other) {
    return x == other.x && y == other.y && z == other.z;
  }

}

class Mesh {
  ArrayList<Integer> V;
  ArrayList<Vertex> G, faceNorms;
  int[] O;
  ArrayList<int[]> colors;
  Vertex[] vertexNorms;
  int num_vertices, num_faces;
  int currCorner = 0;

  Mesh(int v, int f) {
    V = new ArrayList<Integer>();
    G = new ArrayList<Vertex>();
    O = new int[3 * f];
    faceNorms = new ArrayList<Vertex>();
    vertexNorms = new Vertex[v];
    colors = new ArrayList<int[]>();
    num_vertices = v; num_faces = f;
    currCorner = 0;
  }
  
  int NextCorner(int corner) { 
  return (corner / 3) * 3 + (corner + 1) % 3; 
  }
  
  int PrevCorner(int corner) { 
    return NextCorner(NextCorner(corner)); 
  }
    
  
  int OppositeCorner(int corner) {
    int faceStart = (corner / 3) * 3;

    int edgeVertex1 = mesh.V.get((faceStart + (corner + 1) % 3) % mesh.V.size());
    int edgeVertex2 = mesh.V.get((faceStart + (corner + 2) % 3) % mesh.V.size());

    // find the face that shares an edge
    for (int i = 0; i < mesh.V.size(); i += 3) {
        if (i == faceStart) continue; // Skip the current face

        // Check if this face shares the edge
        boolean sharesEdgeVertex1 = mesh.V.get(i) == edgeVertex1 || mesh.V.get(i + 1) == edgeVertex1 || mesh.V.get(i + 2) == edgeVertex1;
        boolean sharesEdgeVertex2 = mesh.V.get(i) == edgeVertex2 || mesh.V.get(i + 1) == edgeVertex2 || mesh.V.get(i + 2) == edgeVertex2;

        if (sharesEdgeVertex1 && sharesEdgeVertex2) {
            // Find the corner opposite to this edge in the adjacent face
            for (int j = 0; j < 3; j++) {
                int oppositeCornerIndex = i + j;
                if (mesh.V.get(oppositeCornerIndex) != edgeVertex1 && mesh.V.get(oppositeCornerIndex) != edgeVertex2) {
                    return oppositeCornerIndex;
                }
            }
        }
    }

    return -1;
  }
  int SwingCorner(int corner) { 
    return NextCorner(O[NextCorner(corner)]); 
  }
  
  // Opposite Table
  void OTable() {
    for (int a = 0; a < V.size(); a++) {
      for (int b = 0; b < V.size(); b++) {
        if (V.get(NextCorner(a)).equals(V.get(PrevCorner(b))) && V.get(PrevCorner(a)).equals(V.get(NextCorner(b)))) {
          O[a] = b;
          O[b] = a;
        }
      }
    }
  }
  
Mesh subdivide(Mesh oldMesh) {
    ArrayList<ArrayList<Integer>> adjVerts = CalcAdjVerts(oldMesh);
    ArrayList<Vertex> G_New = computeEvenVertices(oldMesh, adjVerts);
    ArrayList<Integer> V_New = computeOddVertices(oldMesh, G_New);
    return createNewMesh(G_New, V_New);
}

// Helper Methods

ArrayList<Vertex> computeEvenVertices(Mesh oldMesh, ArrayList<ArrayList<Integer>> adjacentVertices) {
    ArrayList<Vertex> G_New = new ArrayList<Vertex>();
    for (int i = 0; i < oldMesh.G.size(); i++) {
        ArrayList<Integer> adjVerts = adjacentVertices.get(i);
        Vertex NeighborSum = new Vertex(0.0, 0.0, 0.0);
        for (int adjVert : adjVerts) {
            NeighborSum = NeighborSum.add(oldMesh.G.get(adjVert));
        }
        float beta = adjVerts.size() == 3 ? 3.0 / 16.0 : 3.0 / (8.0 * adjVerts.size());
        G_New.add(oldMesh.G.get(i).mult(1 - adjVerts.size() * beta).add(NeighborSum.mult(beta)));
    }
    return G_New;
}


ArrayList<Integer> computeOddVertices(Mesh oldMesh, ArrayList<Vertex> G_New) {
    ArrayList<Integer> V_New = new ArrayList<Integer>();
    ArrayList<int[]> visited = new ArrayList<int[]>(); 
    int[] triangleEdges = new int[3]; //

    int first = 0;
    int second = 0;
    
    for (int j = 0; j < oldMesh.V.size(); j++) {
      int currentVertex = oldMesh.V.get(j);
      int nextVert = oldMesh.V.get(mesh.NextCorner(j));
      int prevVertex = oldMesh.V.get(mesh.PrevCorner(j));
  
      // Define adjacent edges
      int[] firstEdge = {currentVertex, nextVert};
      Arrays.sort(firstEdge);
  
      int[] secondEdge = {currentVertex, prevVertex};
      Arrays.sort(secondEdge);
  
      // Check if these edges have been processed
      boolean ContainsFirst = false;
      boolean ContainsSecond = false;
      for (int k = 0; k < visited.size(); k++) {
        int[] n = visited.get(k);
        if (n[0] == firstEdge[0] && n[1] == firstEdge[1]) ContainsFirst = true;
        if (n[0] == secondEdge[0] && n[1] == secondEdge[1]) ContainsSecond = true;
      }
  
      Vertex lVert = oldMesh.G.get(prevVertex);
      Vertex rVert = oldMesh.G.get(oldMesh.V.get(oldMesh.O[mesh.PrevCorner(j)]));
      
      Vertex currVert = oldMesh.G.get(currentVertex);
      Vertex preVert = oldMesh.G.get(prevVertex);
      Vertex NewVert = currVert.add(oldMesh.G.get(nextVert)).mult(3.0 / 8.0).add(lVert.add(rVert).mult(1.0 / 8.0));
      
      // get index or add
      if (!ContainsFirst) {
        visited.add(firstEdge);
        first = G_New.size();
        G_New.add(NewVert);
      } else {
        for (int k = 0; k < G_New.size(); k++) {
          Vertex check = G_New.get(k);
          if (check.equals(NewVert)) {
            first = k;
            break;
          }
        }
      }
      
      // second edge
      lVert = oldMesh.G.get(nextVert);
      rVert = oldMesh.G.get(oldMesh.V.get(oldMesh.O[mesh.NextCorner(j)]));
      NewVert = currVert.add(preVert).mult(3.0 / 8.0).add(lVert.add(rVert).mult(1.0 / 8.0));
  
      if (!ContainsSecond){
        visited.add(secondEdge);
        second = G_New.size();
        G_New.add(NewVert);
      } else {      
        for (int k = 0; k < G_New.size(); k++) {
          Vertex check = G_New.get(k);
          if (check.equals(NewVert)) {
            second = k;
            break;
          }
        }
      }

      V_New.add(currentVertex);
      V_New.add(first);
      V_New.add(second);
      
      // Handle inner triangle 
      if(j % 3 == 0) triangleEdges[0] = first;
      else if (j % 3 == 1) triangleEdges[1] = first;
      else {
        V_New.add(triangleEdges[0]);
        V_New.add(triangleEdges[1]);
        V_New.add(first);
      }
    }

    return V_New;
}


Mesh createNewMesh(ArrayList<Vertex> newG, ArrayList<Integer> newV) {
    Mesh newMesh = new Mesh(newG.size(), newV.size() / 3);
    newMesh.V = newV;
    newMesh.G = newG;

    newMesh.OTable();
    newMesh.CalcFaceNorms();
    newMesh.CalcVertexNorms();

    return newMesh;
}

ArrayList<ArrayList<Integer>> CalcAdjVerts(Mesh mesh) {
  int numVertices = mesh.G.size();
  ArrayList<ArrayList<Integer>> adjacentVerticesLists = new ArrayList<ArrayList<Integer>>(numVertices);
  boolean[] isVertexProcessed = new boolean[numVertices];

  for (int i = 0; i < numVertices; i++) {
    adjacentVerticesLists.add(new ArrayList<Integer>());
  }

  for (int i = 0; i < mesh.V.size(); i++) {
    int currentVertex = mesh.V.get(i);
    
    if (isVertexProcessed[currentVertex]) continue;

    ArrayList<Integer> adjacentVertices = new ArrayList<Integer>();
    int swingVertex = i;

    // Find all adjacent vertices by swinging around the current vertex
    do {
      int nextVertex = mesh.V.get(mesh.NextCorner(swingVertex));
      adjacentVertices.add(nextVertex);
      swingVertex = mesh.SwingCorner(swingVertex);
    } while (swingVertex != i);

    isVertexProcessed[currentVertex] = true;
    adjacentVerticesLists.set(currentVertex, adjacentVertices);
  }
  return adjacentVerticesLists;
}

  
void CalcFaceNorms() { 
    faceNorms.clear(); // Clear existing normals
    for (int i = 0; i < V.size(); i += 3) {
        // Ensure there are enough vertices for a face
        if (i + 2 < V.size()) {
            Vertex vertexA = G.get(V.get(i));
            Vertex vertexB = G.get(V.get(i + 1));
            Vertex vertexC = G.get(V.get(i + 2));
            Vertex normal = CalcFaceNorm(vertexA, vertexB, vertexC);
            faceNorms.add(normal);
        }
    }
}

Vertex CalcFaceNorm(Vertex a, Vertex b, Vertex c) {
    Vertex edge1 = b.sub(a);
    Vertex edge2 = c.sub(a);
    Vertex norm = edge1.cross(edge2);

    float len = sqrt(norm.x * norm.x + norm.y * norm.y + norm.z * norm.z);

    // Normalize 
    if (len != 0) {
        norm.x /= len;
        norm.y /= len;
        norm.z /= len;
    }
    return norm;
}
    
void CalcVertexNorms() {
    boolean[] normalized = new boolean[G.size()];

    for (int i = 0; i < V.size(); i++) {
      if (normalized[V.get(i)] == false) {
        Vertex norm = faceNorms.get(i / 3);
        int swingV = SwingCorner(i);
        while(swingV != i) {
          norm = norm.add(faceNorms.get(swingV / 3));
          swingV = SwingCorner(swingV);
        }
        
        vertexNorms[V.get(i)] = norm.norm();
        normalized[V.get(i)] = true;
      }
    }
  }
}

// remember old mouse position
void mousePressed()
{
  mouseX_old = mouseX;
  mouseY_old = mouseY;
}

// modify rotation matrix when mouse is dragged
void mouseDragged()
{
  if (!mousePressed)
    return;
  
  float dx = mouseX - mouseX_old;
  float dy = mouseY - mouseY_old;
  dy *= -1;

  float len = sqrt (dx*dx + dy*dy);
  if (len == 0)
      len = 1;
  
  dx /= len;
  dy /= len;
  PMatrix3D rmat = new PMatrix3D();
  rmat.rotate (len * 0.005, dy, dx, 0);
  rot_mat.preApply (rmat);

  mouseX_old = mouseX;
  mouseY_old = mouseY;
}
