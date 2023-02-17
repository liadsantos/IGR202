#ifndef MESH_H
#define MESH_H

#include <glad/glad.h>
#include <vector>
#include <memory>
#include <math.h>
#include <cmath>
#include <cstdlib>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include <map>
#include <set>
#include <random>
#include <iterator>

#include <iostream>

class Mesh {
public:
  virtual ~Mesh();

  struct Edge {
    unsigned int a , b;
    Edge( unsigned int c , unsigned int d ) : a( std::min<unsigned int>(c,d) ) , b( std::max<unsigned int>(c,d) ) {}
    bool operator < ( Edge const & o ) const {   return a < o.a  ||  (a == o.a && b < o.b);  }
    bool operator == ( Edge const & o ) const {   return a == o.a  &&  b == o.b;  }
    
    float calculateDistance(std::vector<glm::vec3> const vecPos) const {
      glm::vec3 vector;

      vector.x = pow(vecPos[b].x - vecPos[a].x, 2);
      vector.y = pow(vecPos[b].y - vecPos[a].y, 2);
      vector.z = pow(vecPos[b].z - vecPos[a].z, 2);

      return sqrt(vector.x + vector.y + vector.z);
    }
    
    glm::vec3 crossProduct(
      std::vector<glm::vec3> const vecPos,
      Edge const &edge1
    ) const {
      glm::vec3 vector1, vector2;
      
      vector1.x = vecPos[a].x - vecPos[b].x;
      vector1.y = vecPos[a].y - vecPos[b].y;
      vector1.z = vecPos[a].z - vecPos[b].z;

      vector2.x = vecPos[edge1.a].x - vecPos[edge1.b].x;
      vector2.y = vecPos[edge1.a].y - vecPos[edge1.b].y;
      vector2.z = vecPos[edge1.a].z - vecPos[edge1.b].z;

      return glm::cross(vector1, vector2);
    }
  };

  const std::vector<glm::vec3> &vertexPositions() const { return _vertexPositions; }
  std::vector<glm::vec3> &vertexPositions() { return _vertexPositions; }

  const std::vector<glm::vec3> &vertexNormals() const { return _vertexNormals; }
  std::vector<glm::vec3> &vertexNormals() { return _vertexNormals; }

  const std::vector<glm::vec2> &vertexTexCoords() const { return _vertexTexCoords; }
  std::vector<glm::vec2> &vertexTexCoords() { return _vertexTexCoords; }

  const std::vector<glm::uvec3> &triangleIndices() const { return _triangleIndices; }
  std::vector<glm::uvec3> &triangleIndices() { return _triangleIndices; }

  /// Compute the parameters of a sphere which bounds the mesh
  void computeBoundingSphere(glm::vec3 &center, float &radius) const;

  void recomputePerVertexNormals(bool angleBased = false);
  void recomputePerVertexTextureCoordinates( );

  void init();
  void initOldGL();
  void render();
  void clear();

  void addPlan(float square_half_side = 1.0f);

  /*float calculateAlpha(unsigned int numTri) {
    float alpha;

    if (numTri > 3) {
      alpha = 0.125f
    }
    else {
      alpha = 0.1875f;
    }
  }*/

  void subdivideLinear() {
    std::vector<glm::vec3> newVertices = _vertexPositions;
    std::vector<glm::uvec3> newTriangles;

    struct Edge {
      unsigned int a , b;
      Edge( unsigned int c , unsigned int d ) : a( std::min<unsigned int>(c,d) ) , b( std::max<unsigned int>(c,d) ) {}
      bool operator < ( Edge const & o ) const {   return a < o.a  ||  (a == o.a && b < o.b);  }
      bool operator == ( Edge const & o ) const {   return a == o.a  &&  b == o.b;  }
    };
    std::map< Edge , unsigned int > newVertexOnEdge;
    for(unsigned int tIt = 0 ; tIt < _triangleIndices.size() ; ++tIt) { // compute each triangle face
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];

      Edge Eab(a,b);
      unsigned int oddVertexOnEdgeEab = 0;
      if( newVertexOnEdge.find( Eab ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ a ] + _vertexPositions[ b ]) / 2.f );
        oddVertexOnEdgeEab = newVertices.size() - 1;
        newVertexOnEdge[Eab] = oddVertexOnEdgeEab;
      }
      else { oddVertexOnEdgeEab = newVertexOnEdge[Eab]; }

      Edge Ebc(b,c);
      unsigned int oddVertexOnEdgeEbc = 0;
      if( newVertexOnEdge.find( Ebc ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ b ] + _vertexPositions[ c ]) / 2.f );
        oddVertexOnEdgeEbc = newVertices.size() - 1;
        newVertexOnEdge[Ebc] = oddVertexOnEdgeEbc;
      }
      else { oddVertexOnEdgeEbc = newVertexOnEdge[Ebc]; }

      Edge Eca(c,a);
      unsigned int oddVertexOnEdgeEca = 0;
      if( newVertexOnEdge.find( Eca ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ c ] + _vertexPositions[ a ]) / 2.f );
        oddVertexOnEdgeEca = newVertices.size() - 1;
        newVertexOnEdge[Eca] = oddVertexOnEdgeEca;
      }
      else { oddVertexOnEdgeEca = newVertexOnEdge[Eca]; }

      // set new triangles :
      newTriangles.push_back( glm::uvec3( a , oddVertexOnEdgeEab , oddVertexOnEdgeEca ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEab , b , oddVertexOnEdgeEbc ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEca , oddVertexOnEdgeEbc , c ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEab , oddVertexOnEdgeEbc , oddVertexOnEdgeEca ) );
    }

    // after that:
    _triangleIndices = newTriangles;
    _vertexPositions = newVertices;
    recomputePerVertexNormals( );
    recomputePerVertexTextureCoordinates( );
  }

  void subdivideLoop2() {
    // Declare new vertices and new triangles. Initialize the new positions for the even vertices with (0,0,0):
    //std::vector<glm::vec3> newVertices( _vertexPositions.size(), glm::vec3(0,0,0) );
    std::vector<glm::vec3> newVertices;
    std::vector<glm::uvec3> newTriangles;

    std::map< Edge , unsigned int > newVertexOnEdge; // this will be useful to find out whether we already inserted an odd vertex or not
    std::map< Edge , std::set< unsigned int > > trianglesOnEdge; // this will be useful to find out if an edge is boundary or not
    std::map< unsigned int, unsigned int > oddVertices;   // vertices in faces that share edges
    std::vector< std::set< unsigned int > > neighboringVertices( _vertexPositions.size() ); // this will be used to store the adjacent vertices, i.e., neighboringVertices[i] will be the list of vertices that are adjacent to vertex i.
    std::vector< bool > evenVertexIsBoundary( _vertexPositions.size() , false );
    std::vector< bool > evenVertexIsShared( _vertexPositions.size() , false );

    // Compute the valences of the even vertices, the neighboring vertices required to update the position of the even vertices, and the boundaries:
    for(unsigned int tIt = 0 ; tIt < _triangleIndices.size() ; ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];

      //Remember the faces shared by the edge
      Edge Eab(a,b);
      Edge Ebc(b,c);
      Edge Eca(c,a);

      if (trianglesOnEdge.find(Eab) == trianglesOnEdge.end()) {
        evenVertexIsShared[a] = true;
        evenVertexIsShared[b] = true;
      } else {
        evenVertexIsShared[a] = false;
        evenVertexIsShared[b] = false;
      }

      if (trianglesOnEdge.find(Ebc) == trianglesOnEdge.end()) {
        evenVertexIsShared[b] = true;
        evenVertexIsShared[c] = true;
      } else {
        evenVertexIsShared[b] = false;
        evenVertexIsShared[c] = false;
      }

      if (trianglesOnEdge.find(Eca) == trianglesOnEdge.end()) {
        evenVertexIsShared[c] = true;
        evenVertexIsShared[a] = true;
      } else {
        evenVertexIsShared[c] = false;
        evenVertexIsShared[a] = false;
      }

      // stores the indices in a set (container)
      neighboringVertices[ a ].insert( b );
      neighboringVertices[ a ].insert( c );
      neighboringVertices[ b ].insert( a );
      neighboringVertices[ b ].insert( c );
      neighboringVertices[ c ].insert( a );
      neighboringVertices[ c ].insert( b );
    }

    // The valence of a vertex is the number of adjacent vertices:
    std::vector< unsigned int > evenVertexValence( _vertexPositions.size() ,0 );
    for( unsigned int v = 0 ; v < _vertexPositions.size() ; ++v ) {
      evenVertexValence[ v ] = neighboringVertices[ v ].size();
    }

    // Set variables to create the new vertices
    float alpha;
    glm::vec3 sum = glm::vec3(0.0,0.0,0.0);
    unsigned int numTriangles;

    // loop for updating the even vertices
    for(unsigned int v = 0; v < _vertexPositions.size(); ++v) {

      // set the number of neighboring triangles
      numTriangles = evenVertexValence[v];

      // Calculate the weights
      if ( numTriangles > 3 )       { alpha = 0.375f / (float)numTriangles; }
      else if ( numTriangles == 3 ) { alpha = 0.1875f; }

      for (std::set<unsigned int> :: iterator it = neighboringVertices[v].begin(); it != neighboringVertices[v].end(); ++it) {
        sum += alpha * _vertexPositions[ *it ];                          // neighboring vertices
      }

      sum += (1 - (evenVertexValence[v] * alpha)) * _vertexPositions[v]; // centroid vertex
      newVertices.push_back(sum);
      sum = glm::vec3(0.0,0.0,0.0);
    }

    // Compute the odd vertices:
    for(unsigned int tIt = 0 ; tIt < _triangleIndices.size() ; ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];

      Edge Eab(a,b);
      unsigned int oddVertexOnEdgeEab = 0;

      // there is no Eab edge, so create vertice
      if( newVertexOnEdge.find( Eab ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ a ] + _vertexPositions[ b ]) / 2.f );
        oddVertexOnEdgeEab = newVertices.size() - 1; // indice of vertex created on the edge
        newVertexOnEdge[Eab] = oddVertexOnEdgeEab;
        oddVertices[oddVertexOnEdgeEab] = c;         // indice of third vertex created from first face of the edge
      }

      // there is Eab edge, so update odd vertice
      else {
        oddVertexOnEdgeEab = newVertexOnEdge[Eab];
        newVertices[oddVertexOnEdgeEab] = newVertices[oddVertexOnEdgeEab] * 3.f / 4.f +             // 3/8*a + 3/8*b
                                          _vertexPositions[c] / 8.f +                               // 1/8*c
                                          _vertexPositions[oddVertices[oddVertexOnEdgeEab]] / 8.f;  // 1/8*d (c from the created Eab face)
        trianglesOnEdge[Eab].insert(a);
        trianglesOnEdge[Eab].insert(b);
      }

      Edge Ebc(b,c);
      unsigned int oddVertexOnEdgeEbc = 0;
      if( newVertexOnEdge.find( Ebc ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ b ] + _vertexPositions[ c ]) / 2.f );
        oddVertexOnEdgeEbc = newVertices.size() - 1;
        newVertexOnEdge[Ebc] = oddVertexOnEdgeEbc;
        oddVertices[oddVertexOnEdgeEbc] = a;
      }
      else {
        oddVertexOnEdgeEbc = newVertexOnEdge[Ebc];
        newVertices[oddVertexOnEdgeEbc] = newVertices[oddVertexOnEdgeEbc] * 3.f / 4.f +
                                          _vertexPositions[a] / 8.f +
                                          _vertexPositions[oddVertices[oddVertexOnEdgeEbc]] / 8.f;
      }

      Edge Eca(c,a);
      unsigned int oddVertexOnEdgeEca = 0;
      if( newVertexOnEdge.find( Eca ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ c ] + _vertexPositions[ a ]) / 2.f );
        oddVertexOnEdgeEca = newVertices.size() - 1;
        newVertexOnEdge[Eca] = oddVertexOnEdgeEca;
        oddVertices[oddVertexOnEdgeEca] = b;
      }
      else {
        oddVertexOnEdgeEca = newVertexOnEdge[Eca];
        newVertices[oddVertexOnEdgeEca] = newVertices[oddVertexOnEdgeEca] * 3.f / 4.f +
                                          _vertexPositions[b] / 8.f +
                                          _vertexPositions[oddVertices[oddVertexOnEdgeEca]] / 8.f;
      }

      // set new triangles :
      newTriangles.push_back( glm::uvec3( a , oddVertexOnEdgeEab , oddVertexOnEdgeEca ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEab , b , oddVertexOnEdgeEbc ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEca , oddVertexOnEdgeEbc , c ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEab , oddVertexOnEdgeEbc , oddVertexOnEdgeEca ) );
    }

    // after that:
    _triangleIndices = newTriangles;
    _vertexPositions = newVertices;
    recomputePerVertexNormals( );
    recomputePerVertexTextureCoordinates( );
  }

  float getAverageEdgeLength() {
    float averageEdgeLength = 0.f;
    std::vector<Edge> edges;

    edges.clear();

    for(unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
      // creates the edges for each triangle on the mesh
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];

      Edge Eab(a,b);
      Edge Ebc(b,c);
      Edge Eca(a,c);

      // add the distance of each edge
      averageEdgeLength += Eab.calculateDistance(_vertexPositions);
      averageEdgeLength += Ebc.calculateDistance(_vertexPositions);
      averageEdgeLength += Eca.calculateDistance(_vertexPositions);

      edges.push_back(Eab);
      edges.push_back(Ebc);
      edges.push_back(Eca);
    }

    return averageEdgeLength /= edges.size();
  }

  std::vector<float> getFaceArea() {
    std::vector<float> area(_triangleIndices.size(), 0.0);
    float s, distA, distB, distC;

    for(unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];

      Edge Eab(a,b);
      Edge Ebc(b,c);
      Edge Eca(c,a);

      s = 0.0, distA = 0.0, distB = 0.0, distC = 0.0;
      distA = Eab.calculateDistance(_vertexPositions);
      distB = Ebc.calculateDistance(_vertexPositions);
      distC = Eca.calculateDistance(_vertexPositions);

      s = 0.5 * (distA + distB + distC);

      area[tIt] = sqrt(s * (s-distA) * (s-distB) * (s-distC));
    }

    return area;
  }

  void getCentroid(std::vector<glm::vec3> &faceCentroid) {
    faceCentroid.clear();                         // initialization
    float aCent = 0.f, bCent = 0.f, cCent = 0.f;  // coordinates for centroid

    for(unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];

      float aCent = (_vertexPositions[a].x + _vertexPositions[b].x + _vertexPositions[c].x) / 3;
      float bCent = (_vertexPositions[a].y + _vertexPositions[b].y + _vertexPositions[c].y) / 3;
      float cCent = (_vertexPositions[a].z + _vertexPositions[b].z + _vertexPositions[c].z) / 3;

      faceCentroid.push_back(glm::vec3(aCent, bCent, cCent));
    }
  }

  std::vector<glm::vec3> getCentroidsNeighboringFaces(
    std::vector<glm::vec3> &faceCentroid,
    std::vector<unsigned int> &indicesOfTriangles,
    glm::vec3 point,
    float sigma_f
  ) {
    float radius = 2.0 * sigma_f;                     // radius of neighbors
    glm::vec3 faceCent;                               // centroid of the considered face in the loop
    std::vector<glm::vec3> neighborsOfPoint;          // storage the neighbors

    neighborsOfPoint.clear();                         // initialization
    for(unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
      faceCent = faceCentroid[tIt];
      float distance = glm::distance(point, faceCent);

      if (distance <= radius && distance > 0.0) {
        neighborsOfPoint.push_back(faceCent);
        indicesOfTriangles.push_back(tIt);
      }
    }

    return neighborsOfPoint;
  }

  glm::vec3 calculateTriangleNormal(Edge &edge1, Edge &edge2) {
    return edge1.crossProduct(_vertexPositions, edge2);
  }

  void mollification(
    std::vector<glm::vec3> &mollifiedNormals,
    std::vector<glm::vec3> &faceCentroid,
    std::vector<float> &faceArea,
    float sigma_f
  ) {
    float sigma_fm = sigma_f / 2.f;                                     // spatial weigth for mollification
    std::vector<std::vector<glm::vec3>> vertNeighFaces;                 // keep the neighboring vertices of all points
    std::vector<unsigned int> indicesOfTriangles;                       // storage the indices of triangles

    for (auto &vert : vertNeighFaces)                                   // initialization of vertices
      vert.clear();

    for (unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];

      Edge Eab(a,b);
      Edge Ebc(b,c);

      glm::vec3 normal = calculateTriangleNormal(Eab, Ebc);             
      glm::vec3 newNormal = glm::vec3(0.0, 0.0, 0.0);
      glm::vec3 centroidFace = faceCentroid[tIt];

      float k = 0.0;                                                    // normalization factor

      indicesOfTriangles.clear();

      vertNeighFaces.push_back(                                         // define the centroids of neighboring faces
        getCentroidsNeighboringFaces(faceCentroid, indicesOfTriangles, centroidFace, sigma_f)
      );

      // Pass of non-robust smoothing using equation (3) without influence weigth
      for (unsigned int i = 0; i < vertNeighFaces[tIt].size(); ++i) {
        unsigned int idx = indicesOfTriangles[i];

        float dis = glm::distance(centroidFace, faceCentroid[idx]);
        float f = std::exp(-dis * dis / (2 * sigma_fm * sigma_fm));     // spatial weight
        float aq = faceArea[idx];

        newNormal += faceCentroid[idx] * aq * f;
        k += aq * f;
      }
      newNormal /= k;

      mollifiedNormals.push_back(newNormal);                            // push the new mollified normals
    }
  }

  void calculateTriangleNormals(std::vector<glm::vec3> &mollifiedNormals) {
    for(unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];

      Edge Eab(a,b);
      Edge Ebc(b,c);

      glm::vec3 n = calculateTriangleNormal(Eab, Ebc);

      mollifiedNormals.push_back(calculateTriangleNormal(Eab, Ebc));
    }
  }
    
  void addNoise(){
    for(unsigned int i = 0 ; i < _vertexPositions.size(); ++i) {
      _vertexPositions[i].x += ((rand() % 100) - 50) * 0.001;
      _vertexPositions[i].y += ((rand() % 100) - 50) * 0.001;
      _vertexPositions[i].z += ((rand() % 100) - 50) * 0.001;
    }
    
    recomputePerVertexNormals( );
    recomputePerVertexTextureCoordinates( );
  }

  /* Estimate the new vertex positions
   * receives the following parameters:
   * sigma_fp: value of sigma_f from table 1
   * sigma_gp: value of sigma_g from table 1
   * position[0,1,2]: for different meshes, different sigmas can be tested. The user specify which one he wants here
   */
  void robustEstimation(std::vector<float> sigma_fp, std::vector<float> sigma_gp, unsigned int position, bool mollify) {
    // ----------- Useful variables
    float sigma_f;                            // spatial weight gaussian
    float sigma_g;                            // influence weight gaussian
    float k;                                  // normalization factor
    float f, g;                               // spatial and influence weights
    float meanEdgeLength;                     // ||e||
    std::vector<float> faceArea;              // area of the face
    std::vector<glm::vec3> faceCentroid;      // centroid of the face
    std::vector<glm::vec3> newVertices;       // push the changed vertices
    std::vector<glm::vec3> mollifiedNormals;  // smoothed normals

    // ----------- Smoothing and denoising
    meanEdgeLength = getAverageEdgeLength();  // take a cube
    sigma_f = sigma_fp[position] * meanEdgeLength;
    sigma_g = sigma_gp[position] * meanEdgeLength;
    faceArea = getFaceArea();
    getCentroid(faceCentroid);

    // std::cout << "sigma f: " << sigma_f << "\t\t";
    // std::cout << "sigma f_ant: " << sigma_fp[position] << "\t\t";
    // std::cout << "mean edge: " << meanEdgeLength << "\t\t";

    // mollification
    mollifiedNormals.clear();
    if (mollify)
      mollification(mollifiedNormals, faceCentroid, faceArea, sigma_f);
    else 
      calculateTriangleNormals(mollifiedNormals);

    // robust estimation of vertices
    newVertices.clear();
    std::cout << "Point: " << _vertexPositions[11].x << "\t\t" << _vertexPositions[11].y << "\t\t" << _vertexPositions[11].z << "\n";

    for (unsigned int vIt = 0; vIt < _vertexPositions.size(); ++vIt) {
      // point to be estimated
      glm::vec3 p = _vertexPositions[vIt]; 

      // indices of each neighbor triangle
      std::vector<unsigned int> indicesOfTriangles;                             
      indicesOfTriangles.clear();

      // get the neighboring faces
      std::vector<glm::vec3> vertNeighFaces = getCentroidsNeighboringFaces(
        faceCentroid, indicesOfTriangles, p, sigma_f
      );

      // parameters for equation (3) and (4)
      glm::vec3 new_p = glm::vec3(0.0, 0.0, 0.0);                                          
      float k = 0.0;                                                    

      for (unsigned int i = 0; i < vertNeighFaces.size(); ++i) {
        unsigned int idx = indicesOfTriangles[i];
        
        glm::vec3 normal = mollifiedNormals[idx];                                                       // normal of the neighbor face
        glm::vec3 cent = faceCentroid[idx];                                                             // centroid of the neighbor face
        glm::vec3 predictor_g = p - glm::dot((p-cent), normal) * normal;                                // PI_q(p)                                           

        float dist_spacial = glm::distance(p, faceCentroid[idx]);                                       // ||p - c_q||
        float weight_spacial = std::exp(-dist_spacial * dist_spacial / (2 * sigma_f * sigma_f));        // f

        float dist_influence = glm::distance(p, predictor_g);                                           // ||p - PI_q(p)||
        float weight_influence = std::exp(-dist_influence * dist_influence / (2 * sigma_g * sigma_g));  // g

        float area = faceArea[idx];

        new_p += predictor_g * area * weight_spacial * weight_influence;
        k += area * weight_spacial * weight_influence;
      }
      new_p /= k;
      newVertices.push_back(new_p);   
      std::cout << "k: " << k << "\n";
      
    }
    
    std::cout << "Point: " << newVertices[11].x << "\t\t" << newVertices[11].y << "\t\t" << newVertices[11].z << "\n";

    _vertexPositions = newVertices;
    std::cout << "Finished!\n";
    
    recomputePerVertexNormals( );
    recomputePerVertexTextureCoordinates( );
  }

  void subdivideLoop() {
    //subdivideLinear();
    subdivideLoop2();
  }

  void applyMeshDenoising(
    std::vector<float> sigma_fp, 
    std::vector<float> sigma_gp, 
    unsigned int position,
    bool mollify
  ) {
    robustEstimation(sigma_fp, sigma_gp, position, mollify);
  }

private:
  std::vector<glm::vec3> _vertexPositions;
  std::vector<glm::vec3> _vertexNormals;
  std::vector<glm::vec2> _vertexTexCoords;
  std::vector<glm::uvec3> _triangleIndices;

  GLuint _vao = 0;
  GLuint _posVbo = 0;
  GLuint _normalVbo = 0;
  GLuint _texCoordVbo = 0;
  GLuint _ibo = 0;
};

// utility: loader
void loadOFF(const std::string &filename, std::shared_ptr<Mesh> meshPtr);

#endif  // MESH_H
