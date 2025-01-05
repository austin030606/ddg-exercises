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
    MeshSubset star_subset;
    star_subset.addVertices(subset.vertices);
    star_subset.addEdges(subset.edges);
    star_subset.addFaces(subset.faces);
    Vector<size_t> star_edges_vec, 
                   star_faces_vec, 
                   subset_vertex_vec = buildVertexVector(subset), 
                   subset_edge_vec = buildEdgeVector(subset),
                   subset_face_vec = buildFaceVector(subset);
    
    star_edges_vec = A0 * subset_vertex_vec + subset_edge_vec;
    for (size_t i = 0; i < star_edges_vec.size(); i++) {
        if (star_edges_vec[i] != 0) {
            star_edges_vec[i] = 1;
            star_subset.addEdge(i);
        }
    }

    star_faces_vec = A1 * star_edges_vec + subset_face_vec;
    for (size_t i = 0; i < star_faces_vec.size(); i++) {
        if (star_faces_vec[i] != 0) {
            star_subset.addFace(i);
        }
    }
    return star_subset;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    MeshSubset closure_subset;
    closure_subset.addVertices(subset.vertices);
    closure_subset.addEdges(subset.edges);
    closure_subset.addFaces(subset.faces);

    Vector<size_t> closure_edges_vec, 
                   closure_vertices_vec, 
                   subset_vertex_vec = buildVertexVector(subset), 
                   subset_edge_vec = buildEdgeVector(subset),
                   subset_face_vec = buildFaceVector(subset);

    closure_edges_vec = A1.transpose() * subset_face_vec + subset_edge_vec;
    for (size_t i = 0; i < closure_edges_vec.size(); i++) {
        if (closure_edges_vec[i] != 0) {
            closure_edges_vec[i] = 1;
            closure_subset.addEdge(i);
        }
    }

    closure_vertices_vec = A0.transpose() * closure_edges_vec + subset_vertex_vec;
    for (size_t i = 0; i < closure_vertices_vec.size(); i++) {
        if (closure_vertices_vec[i] != 0) {
            closure_subset.addVertex(i);
        }
    }
    return closure_subset;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    MeshSubset closure_of_star = closure(star(subset)),
               star_of_closure = star(closure(subset));

    closure_of_star.deleteSubset(star_of_closure);
    return closure_of_star;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    return subset.equals(closure(subset)); // placeholder
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
    if (!isComplex(subset))
        return -1;
    
    if (!subset.faces.empty()) {
        MeshSubset faces_subset, closure_subset;
        faces_subset.addFaces(subset.faces);
        closure_subset = closure(faces_subset);
        if (closure_subset.equals(subset)) {
            return 2;
        }
    } else if (!subset.edges.empty()) {
        MeshSubset edges_subset, closure_subset;
        edges_subset.addEdges(subset.edges);
        closure_subset = closure(edges_subset);
        if (closure_subset.equals(subset)) {
            return 1;
        }
    } else if (!subset.vertices.empty()) {
        MeshSubset vertices_subset, closure_subset;
        vertices_subset.addVertices(subset.vertices);
        closure_subset = closure(vertices_subset);
        if (closure_subset.equals(subset)) {
            return 0;
        }
    }
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
    if (isPureComplex(subset) == -1) {
        std::cout << "selected subset is not pure\n";
        return subset;
    }
    MeshSubset boundaryFaces; // the set of all simplices that are proper faces of exactly one simplex of subset
    Vector<size_t> nface_per_vertex_vec,
                   nedge_per_vertex_vec,
                   nsimplex_per_vertex_vec,
                   nface_per_edge_vec,
                   subset_vertex_vec = buildVertexVector(subset), 
                   subset_edge_vec = buildEdgeVector(subset),
                   subset_face_vec = buildFaceVector(subset);

    nedge_per_vertex_vec = A0.transpose() * subset_edge_vec;
    
    SparseMatrix<size_t> A1A0 = A1 * A0;
    for (size_t i = 0; i < A1A0.nonZeros(); i++) {
        A1A0.valuePtr()[i] = 1;
    }

    nface_per_vertex_vec = A1A0.transpose() * subset_face_vec;
    nsimplex_per_vertex_vec = nedge_per_vertex_vec + nface_per_vertex_vec;
    
    nface_per_edge_vec = A1.transpose() * subset_face_vec;

    for (size_t i = 0; i < nsimplex_per_vertex_vec.size(); i++) {
        if (nsimplex_per_vertex_vec[i] == 1) {
            boundaryFaces.addVertex(i);
        }
    }
    for (size_t i = 0; i < nface_per_edge_vec.size(); i++) {
        if (nface_per_edge_vec[i] == 1) {
            boundaryFaces.addEdge(i);
        }
    }

    return closure(boundaryFaces);
}