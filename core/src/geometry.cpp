// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {

    // TODO
    // assuming triangle mesh
    Halfedge next_he = he.next();
    Vector3 v_a = inputVertexPositions[he.tipVertex()] - inputVertexPositions[next_he.tipVertex()],
            v_b = inputVertexPositions[he.tailVertex()] - inputVertexPositions[next_he.tipVertex()];
    
    Vector3 cross_product = cross(v_a, v_b);

    return dot(v_a, v_b) / cross_product.norm(); 
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {

    // TODO
    double res = 0.0;
    for (Face f: v.adjacentFaces()) {
        res += faceArea(f);
    }
    return res / 3;
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {

    // TODO
    Vertex vertex_c = c.vertex(), vertex_a = c.halfedge().tipVertex(), vertex_b = c.halfedge().next().tipVertex();
    Vector3 v_a = inputVertexPositions[vertex_a] - inputVertexPositions[vertex_c],
            v_b = inputVertexPositions[vertex_b] - inputVertexPositions[vertex_c];
   
    double cos_value = (dot(v_a, v_b)) / (norm(v_a) * norm(v_b));
    if (cos_value > 1.0) {
        cos_value = 1.0;
    } else if (cos_value < -1.0) {
        cos_value = -1.0;
    }

    return acos(cos_value);
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {

    // TODO
    Vertex v_i = he.tailVertex(), v_j = he.tipVertex(), v_k = he.next().tipVertex(), v_l = he.twin().next().tipVertex();
    
    Vector3 e_ij = inputVertexPositions[v_j] - inputVertexPositions[v_i],
            e_ik = inputVertexPositions[v_k] - inputVertexPositions[v_i],
            e_il = inputVertexPositions[v_l] - inputVertexPositions[v_i];

    Vector3 v_ijk = cross(e_ij, e_ik), v_jil = cross(e_il, e_ij);
    Vector3 N_ijk = v_ijk / norm(v_ijk), N_jil = v_jil / norm(v_jil);

    return atan2(dot(e_ij / norm(e_ij), cross(N_ijk, N_jil)), dot(N_ijk, N_jil)); // placeholder
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {

    // TODO
    Vector3 N = Vector3{0., 0., 0.};
    for (Corner c: v.adjacentCorners()) {
        Vertex v_i = c.vertex(), v_j = c.halfedge().tipVertex(), v_k = c.halfedge().next().tipVertex();
        Vector3 e_ij = inputVertexPositions[v_j] - inputVertexPositions[v_i],
                e_ik = inputVertexPositions[v_k] - inputVertexPositions[v_i];
        
        Vector3 normal = cross(e_ij, e_ik);
        normal /= norm(normal);
        N += normal;
    }
    return N / norm(N);
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {

    // TODO
    Vector3 N = Vector3{0., 0., 0.};
    for (Corner c: v.adjacentCorners()) {
        Vertex v_i = c.vertex(), v_j = c.halfedge().tipVertex(), v_k = c.halfedge().next().tipVertex();
        Vector3 e_ij = inputVertexPositions[v_j] - inputVertexPositions[v_i],
                e_ik = inputVertexPositions[v_k] - inputVertexPositions[v_i];
        
        Vector3 normal = cross(e_ij, e_ik);
        normal /= norm(normal);
        N += angle(c) * normal;
    }
    return N / norm(N);
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {

    // TODO
    Vector3 N = Vector3{0., 0., 0.};
    for (Corner c: v.adjacentCorners()) {
        Vertex v_i = c.vertex(), v_j = c.halfedge().tipVertex(), v_k = c.halfedge().next().tipVertex();
        Vector3 e_ij = inputVertexPositions[v_j] - inputVertexPositions[v_i],
                e_ik = inputVertexPositions[v_k] - inputVertexPositions[v_i];
        
        Vector3 cross_product = cross(e_ij, e_ik);
        N += cross_product / (pow(norm(e_ij) * norm(e_ik), 2));
    }
    return N / norm(N);
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {

    // TODO
    Vector3 N = Vector3{0., 0., 0.};
    for (Corner c: v.adjacentCorners()) {
        Vertex v_i = c.vertex(), v_j = c.halfedge().tipVertex(), v_k = c.halfedge().next().tipVertex();
        Vector3 e_ij = inputVertexPositions[v_j] - inputVertexPositions[v_i],
                e_ik = inputVertexPositions[v_k] - inputVertexPositions[v_i];
        
        Vector3 cross_product = cross(e_ij, e_ik);
        Vector3 normal = cross_product / norm(cross_product);
        double face_area = 0.5 * norm(cross_product);
        N += face_area * normal;
    }
    return N / norm(N);
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {

    // TODO
    Vector3 N = Vector3{0., 0., 0.};
    for (Halfedge he: v.outgoingHalfedges()) {
        Vertex v_i = he.tailVertex(), v_j = he.tipVertex();
        Vector3 e_ij = inputVertexPositions[v_j] - inputVertexPositions[v_i];
        
        N += (dihedralAngle(he) / norm(e_ij)) * e_ij;
    }
    return N / norm(N);
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {

    // TODO
    Vector3 N = Vector3{0., 0., 0.};
    for (Halfedge he: v.outgoingHalfedges()) {
        Vertex v_i = he.tailVertex(), v_j = he.tipVertex();
        Vector3 e_ij = inputVertexPositions[v_j] - inputVertexPositions[v_i];
        
        N += (cotan(he) + cotan(he.twin())) * e_ij;
    }
    return N / norm(N);
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {

    // TODO
    double total_interior_angle = 0.0;
    for (Corner c: v.adjacentCorners()) {
        total_interior_angle += angle(c);
    }
    return 2 * PI - total_interior_angle;
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {

    // TODO
    double total_angle_defect = 0.0;
    for (Vertex v: mesh.vertices()) {
        total_angle_defect += angleDefect(v);
    }
    return total_angle_defect;
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {

    // TODO
    double mean_curavture = 0.0;
    for (Halfedge he: v.outgoingHalfedges()) {
        Vertex v_i = he.tailVertex(), v_j = he.tipVertex();
        Vector3 e_ij = inputVertexPositions[v_j] - inputVertexPositions[v_i];
        
        mean_curavture += dihedralAngle(he) * norm(e_ij);
    }
    mean_curavture *= 0.5;
    return mean_curavture; // placeholder
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {

    // TODO
    double dual_area = 0.0;
    for (Halfedge he: v.outgoingHalfedges()) {
        Vertex v_i = he.tailVertex(), v_j = he.tipVertex(), v_k = he.next().tipVertex();
        Vector3 e_ij = inputVertexPositions[v_j] - inputVertexPositions[v_i],
                e_ik = inputVertexPositions[v_k] - inputVertexPositions[v_i];
        
        dual_area += pow(norm(e_ik), 2) * cotan(he.next().next()) + pow(norm(e_ij), 2) * cotan(he);
    }
    dual_area *= 0.125;
    return dual_area;
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {

    // TODO
    double circumcentric_dual_area = circumcentricDualArea(v);
    double K = angleDefect(v) / circumcentric_dual_area, H = scalarMeanCurvature(v) / circumcentric_dual_area;

    double diff = sqrt(H * H - K);
    return std::make_pair(H + diff, H - diff);
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    // TODO
    return identityMatrix<std::complex<double>>(1); // placeholder
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral