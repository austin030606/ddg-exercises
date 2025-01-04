// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();

    std::vector<Eigen::Triplet<size_t>> tripletList;

    for (Edge e : mesh->edges()) {
        size_t edge_index = geometry->edgeIndices[e],
               first_vertex_index = geometry->vertexIndices[e.firstVertex()],
               second_vertex_index = geometry->vertexIndices[e.secondVertex()];

        tripletList.push_back(Eigen::Triplet<size_t>(edge_index, first_vertex_index, 1));
        tripletList.push_back(Eigen::Triplet<size_t>(edge_index, second_vertex_index, 1));
    }

    SparseMatrix<size_t> A0(mesh->nEdges(),mesh->nVertices());
    A0.setFromTriplets(tripletList.begin(), tripletList.end());
    A0.makeCompressed();
    return A0;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    std::vector<Eigen::Triplet<size_t>> tripletList;

    for (Face f : mesh->faces()) {
        size_t face_index = geometry->faceIndices[f];

        for (Edge e : f.adjacentEdges()) {
            size_t edge_index = geometry->edgeIndices[e];
            tripletList.push_back(Eigen::Triplet<size_t>(face_index, edge_index, 1));
        }
    }

    SparseMatrix<size_t> A1(mesh->nFaces(),mesh->nEdges());
    A1.setFromTriplets(tripletList.begin(), tripletList.end());
    A1.makeCompressed();
    return A1;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> vertex_vec;
    vertex_vec.setZero((mesh->nVertices()));
    for (size_t v_index : subset.vertices) {
        vertex_vec[v_index] = 1;
    }
    return vertex_vec;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> edge_vec;
    edge_vec.setZero(mesh->nEdges());
    for (size_t e_index : subset.edges) {
        edge_vec[e_index] = 1;
    }
    return edge_vec;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> face_vec;
    face_vec.setZero(mesh->nFaces());
    for (size_t f_index : subset.faces) {
        face_vec[f_index] = 1;
    }
    return face_vec;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO
    MeshSubset star;
    star.addVertices(subset.vertices);
    star.addEdges(subset.edges);
    star.addFaces(subset.faces);
    Vector<size_t> star_edges_vec, 
                   star_faces_vec, 
                   subset_vertex_vec = buildVertexVector(subset), 
                   subset_edge_vec = buildEdgeVector(subset),
                   subset_face_vec = buildFaceVector(subset);
    
    star_edges_vec = A0 * subset_vertex_vec + subset_edge_vec;
    for (size_t i = 0; i < star_edges_vec.size(); i++) {
        if (star_edges_vec[i] != 0) {
            star_edges_vec[i] = 1;
            star.addEdge(i);
        }
    }

    star_faces_vec = A1 * star_edges_vec + subset_face_vec;
    for (size_t i = 0; i < star_faces_vec.size(); i++) {
        if (star_faces_vec[i] != 0) {
            star.addFace(i);
        }
    }
    return star;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    MeshSubset closure;
    closure.addVertices(subset.vertices);
    closure.addEdges(subset.edges);
    closure.addFaces(subset.faces);
    return closure;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    return false; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}