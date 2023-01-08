#ifndef MESH_H
#define MESH_H

#include <glad/glad.h>
#include <vector>
#include <memory>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include <map>
#include <set>

class Mesh {
public:
  virtual ~Mesh();

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
    for(unsigned int tIt = 0 ; tIt < _triangleIndices.size() ; ++tIt) {
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
    std::vector<glm::vec3> newVertices( _vertexPositions.size(), glm::vec3(0,0,0) );
    std::vector<glm::uvec3> newTriangles;

    struct Edge {
      unsigned int a , b;
      Edge( unsigned int c , unsigned int d ) : a( std::min<unsigned int>(c,d) ) , b( std::max<unsigned int>(c,d) ) {}
      bool operator < ( Edge const & o ) const {   return a < o.a  ||  (a == o.a && b < o.b);  }
      bool operator == ( Edge const & o ) const {   return a == o.a  &&  b == o.b;  }
    };

    std::map< Edge , unsigned int > newVertexOnEdge; // this will be useful to find out whether we already inserted an odd vertex or not
    std::map< Edge , std::set< unsigned int > > trianglesOnEdge; // this will be useful to find out if an edge is boundary or not
    std::vector< std::set< unsigned int > > neighboringVertices( _vertexPositions.size() ); // this will be used to store the adjacent vertices, i.e., neighboringVertices[i] will be the list of vertices that are adjacent to vertex i.
    std::vector< bool > evenVertexIsBoundary( _vertexPositions.size() , false );

    // Compute the valences of the even vertices, the neighboring vertices required to update the position of the even vertices, and the boundaries:
    for(unsigned int tIt = 0 ; tIt < _triangleIndices.size() ; ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];

      //TODO: Remember the faces shared by the edge (DO IT LATER)

      neighboringVertices[ a ].insert( b );
      neighboringVertices[ a ].insert( c );
      neighboringVertices[ b ].insert( a );
      neighboringVertices[ b ].insert( c );
      neighboringVertices[ c ].insert( a );
      neighboringVertices[ c ].insert( b );
    }

    // The valence of a vertex is the number of adjacent vertices:
    std::vector< unsigned int > evenVertexValence( _vertexPositions.size() , 0 );
    for( unsigned int v = 0 ; v < _vertexPositions.size() ; ++v ) {
      evenVertexValence[ v ] = neighboringVertices[ v ].size();
    }

    // Set variables to create the new vertices
    float alpha = 0;
    int sum = 0;
    int numTriangles = _triangleIndices.size();

    // Calculate the weights
    if (numTriangles > 3) { alpha = 3/8 * numTriangles; }
    else if (numTriangles == 3) { alpha = 3/16; }

    // loop for updating the even vertices
    for(unsigned int v = 0; v < _vertexPositions.size(); ++v) {
      for (unsigned int v_neigh = 0; v_neigh < neighboringVertices[v].size(); ++v_neigh) {
        sum += (alpha / numTriangles) * neighboringVertices[v_neigh];    // neighboring vertices
      }
      sum += (1 - alpha) * _vertexPositions[v];  // centroid vertex
      newVertices[v] = sum;   // update the vertex
      sum = 0;
    }

    // Compute the odd vertices:
    for(unsigned int tIt = 0 ; tIt < _triangleIndices.size() ; ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];

      Edge Eab(a,b);
      unsigned int oddVertexOnEdgeEab = 0;
      if( newVertexOnEdge.find( Eab ) == newVertexOnEdge.end() ) {
        newVertices.push_back( glm::vec3(0,0,0) );
        oddVertexOnEdgeEab = newVertices.size() - 1;
        newVertexOnEdge[Eab] = oddVertexOnEdgeEab;
      }
      else { oddVertexOnEdgeEab = newVertexOnEdge[Eab]; }

      //TODO: Update odd vertices


      Edge Ebc(b,c);
      unsigned int oddVertexOnEdgeEbc = 0;
      if( newVertexOnEdge.find( Ebc ) == newVertexOnEdge.end() ) {
        newVertices.push_back( glm::vec3(0,0,0) );
        oddVertexOnEdgeEbc = newVertices.size() - 1;
        newVertexOnEdge[Ebc] = oddVertexOnEdgeEbc;
      }
      else { oddVertexOnEdgeEbc = newVertexOnEdge[Ebc]; }

      //TODO: Update odd vertices


      Edge Eca(c,a);
      unsigned int oddVertexOnEdgeEca = 0;
      if( newVertexOnEdge.find( Eca ) == newVertexOnEdge.end() ) {
        newVertices.push_back( glm::vec3(0,0,0) );
        oddVertexOnEdgeEca = newVertices.size() - 1;
        newVertexOnEdge[Eca] = oddVertexOnEdgeEca;
      }
      else { oddVertexOnEdgeEca = newVertexOnEdge[Eca]; }

      //TODO: Update odd vertices


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

  void subdivideLoop() {
    subdivideLinear();
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
