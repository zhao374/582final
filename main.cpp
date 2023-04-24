#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>

// Function to read the vertex positions from a text file
void read_vertex_positions(const std::string &filename, Eigen::MatrixXd &V) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<Eigen::RowVector3d> vertices;
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        double x, y, z;
        if (!(iss >> x >> y >> z)) {
            break;
        }
        vertices.push_back(Eigen::RowVector3d(x, y, z));
    }
    infile.close();

    V.resize(vertices.size(), 3);
    for (size_t i = 0; i < vertices.size(); ++i) {
        V.row(i) = vertices[i];
    }
}
void updateMaterialFrame(Eigen::MatrixX3d prev_d3, Eigen::MatrixX3d d3)
{
    //TODO
}
void computeGradientAndHessian(Eigen::VectorXd& gradient,
        Eigen::SparseMatrix<double>& hessian,
        Eigen::MatrixX3d& d3,
        Eigen::VectorXd& twist)
{
    //TODO
}
void updateFrameTheta(Eigen::VectorXd& gradient)
{
    //TODO
}
Eigen::Vector3d parallelTransport(Eigen::Vector3d v, Eigen::Vector3d r1, Eigen::Vector3d r2)
{
    //TODO
}
Eigen::VectorXd getTwist(Eigen::MatrixX3d& d2, Eigen::MatrixX3d& d3)
{
    //TODO
}
double applyBendingForce(Eigen::VectorXd& gradient,
        std::vector<Eigen::Triplet<double>>& hessian,
        Eigen::MatrixX3d& d3,
        Eigen::VectorXd& bending_force)
{
    //TODO
}
double applyTwistingForce(Eigen::VectorXd& gradient,
        std::vector<Eigen::Triplet<double>>& hessian,
        Eigen::VectorXd& twist,
        Eigen::VectorXd& twisting_force)
{
    //TODO
}
int main(int argc, char *argv[]) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi E;

    // Read vertex positions from the text file
    read_vertex_positions("/Users/zhaoyanpeng/582final/vertices.txt", V);

    // Create edges between the consecutive vertices
    int num_vertices = V.rows();
    E.resize(num_vertices, 2);
    for (int i = 0; i < num_vertices-1; ++i) {
        E(i, 0) = i;
        E(i, 1) = (i + 1) % num_vertices;
    }

    // Visualize the points and edges using the libigl viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.data().add_points(V,Eigen::RowVector3d(1,0,0));
    for (unsigned i=0;i<E.rows(); ++i)
      viewer.data().add_edges
      (
        V.row(E(i,0)),
        V.row(E(i,1)),
        Eigen::RowVector3d(1,0,0)
      );

    viewer.launch();

    return 0;
}


