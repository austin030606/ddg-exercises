// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripletList;
    for (Vertex v: mesh.vertices()) {
        size_t index = v.getIndex();
        tripletList.push_back(Eigen::Triplet<double>(index, index, barycentricDualArea(v)));
    }

    SparseMatrix<double> hstar_0(mesh.nVertices(), mesh.nVertices());
    hstar_0.setFromTriplets(tripletList.begin(), tripletList.end());
    hstar_0.makeCompressed();
    return hstar_0;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripletList;
    for (Edge e: mesh.edges()) {
        size_t index = e.getIndex();
        Halfedge he = e.halfedge();
        tripletList.push_back(Eigen::Triplet<double>(index, index, (cotan(he) + cotan(he.twin())) / 2));
    }

    SparseMatrix<double> hstar_1(mesh.nEdges(), mesh.nEdges());
    hstar_1.setFromTriplets(tripletList.begin(), tripletList.end());
    hstar_1.makeCompressed();
    return hstar_1;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripletList;
    for (Face f: mesh.faces()) {
        size_t index = f.getIndex();
        tripletList.push_back(Eigen::Triplet<double>(index, index, 1 / faceArea(f)));
    }

    SparseMatrix<double> hstar_2(mesh.nFaces(), mesh.nFaces());
    hstar_2.setFromTriplets(tripletList.begin(), tripletList.end());
    hstar_2.makeCompressed();
    return hstar_2;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

} // namespace surface
} // namespace geometrycentral