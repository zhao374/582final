#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>

struct Rod {
    Eigen::MatrixXd V;
    Eigen::MatrixXi E;
    Eigen::VectorXd L;
};
struct RodProperties {
    double bending_stiffness;
    double twisting_stiffness;
    double max_bending_angle;
    double max_twisting_angle;
    double timestep;
    int num_iterations;
};


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

void create_edges_and_lengths(const Eigen::MatrixXd &V, Eigen::MatrixXi &E, Eigen::VectorXd &lengths) {
    int num_vertices = V.rows();
    E.resize(num_vertices - 1, 2);
    lengths.resize(num_vertices - 1);
    for (int i = 0; i < num_vertices - 1; ++i) {
        E(i, 0) = i;
        E(i, 1) = i + 1;
        lengths(i) = (V.row(E(i, 1)) - V.row(E(i, 0))).norm();
    }
}

void compute_bishop_frames(const Rod &rod, std::vector<Eigen::Matrix3d> &bishop_frames) {
    int num_edges = rod.E.rows();
    bishop_frames.resize(num_edges);
    
    Eigen::Vector3d ref(0, 0, 1);
    Eigen::Vector3d tangent, binormal, normal;
    
    for (int i = 0; i < num_edges; ++i) {
        tangent = (rod.V.row(rod.E(i, 1)) - rod.V.row(rod.E(i, 0))).normalized();
        
        if (i == 0) {
            binormal = tangent.cross(ref).normalized();
        } else {
            Eigen::Vector3d prev_tangent = (rod.V.row(rod.E(i - 1, 1)) - rod.V.row(rod.E(i - 1, 0))).normalized();
            Eigen::Quaterniond rotation = Eigen::Quaterniond::FromTwoVectors(prev_tangent, tangent);
            binormal = rotation * binormal;
        }
        normal = tangent.cross(binormal).normalized();
        bishop_frames[i] << normal, binormal, tangent;
    }
}

void compute_material_frames(const std::vector<Eigen::Matrix3d> &bishop_frames, const Eigen::VectorXd &thetas, std::vector<Eigen::Matrix3d> &material_frames) {
    int num_frames = bishop_frames.size();
    material_frames.resize(num_frames);
    
    for (int i = 0; i < num_frames; ++i) {
        const Eigen::Matrix3d &bishop_frame = bishop_frames[i];
        double theta = thetas(i);
        
        Eigen::Vector3d tangent = bishop_frame.col(0); // The tangent is the first column of the bishop frame
        Eigen::AngleAxisd rotation(theta, tangent); // Create a rotation around the tangent by angle theta
        material_frames[i] = bishop_frame * rotation;
    }
}


void compute_forces(const Rod &rod, const RodProperties &rod_properties, Eigen::MatrixXd &bend,std::vector<double> &twist) {
    
    int num_vertices = rod.V.rows();
    twist.resize(num_vertices);
    bend.setZero(num_vertices, 3);
    Eigen::VectorXd theta;
    theta.setZero(num_vertices-1);
    theta(0) = -M_PI / 2;
    theta(num_vertices-2) = M_PI / 2;
    std::vector<Eigen::Matrix3d> bishop_frame;
    std::vector<Eigen::Matrix3d> material_frame;
    compute_bishop_frames(rod, bishop_frame);
    compute_material_frames(bishop_frame, theta, material_frame);
    for (int i = 1; i < num_vertices - 1; ++i) {
        Eigen::Vector3d edge1 = rod.V.row(rod.E(i - 1, 1)) - rod.V.row(rod.E(i - 1, 0));
        Eigen::Vector3d edge2 = rod.V.row(rod.E(i, 1)) - rod.V.row(rod.E(i, 0));
        Eigen::Vector3d kappa_b_i = 2 * (edge1.cross(edge2)) / (edge1.norm() * edge2.norm() + edge1.dot(edge2));
        Eigen::Vector3d bending_energy = rod_properties.bending_stiffness * kappa_b_i;
        bend.row(i) = bending_energy;
        
        // Compute twisting force
        double twisting_angle = theta(i)-theta(i - 1);
        double twisting_energy= rod_properties.twisting_stiffness * (twisting_angle);
        twist[i]=twisting_energy;
        
        
    }
}

void compute_material_frame_angles(const Rod &rod, Eigen::VectorXd &theta) {
    int num_edges = rod.E.rows();
    theta.resize(num_edges - 1);
    
    // Initialize material frames (dummy example, replace with actual frames)
    std::vector<Eigen::Matrix3d> material_frames(num_edges);
    for (int i = 0; i < num_edges; ++i) {
        material_frames[i].setIdentity();
    }
    
    for (int i = 0; i < num_edges - 1; ++i) {
        Eigen::Vector3d edge1 = rod.V.row(rod.E(i, 1)) - rod.V.row(rod.E(i, 0));
        Eigen::Vector3d edge2 = rod.V.row(rod.E(i + 1, 1)) - rod.V.row(rod.E(i + 1, 0));
        
        Eigen::Vector3d rotation_axis = edge1.cross(edge2).normalized();
        
        // Compute rotation matrix between the consecutive edges
        Eigen::Matrix3d rotation_matrix = Eigen::Quaterniond::FromTwoVectors(edge1, edge2).toRotationMatrix();
        
        // Compute the angle between material frames
        Eigen::Matrix3d relative_rotation = material_frames[i + 1] * rotation_matrix.transpose() * material_frames[i];
        Eigen::AngleAxisd angle_axis(relative_rotation);
        theta(i) = angle_axis.angle() * rotation_axis.dot(angle_axis.axis());
    }
}
/**
 void compute_forces(const Rod &rod, const RodProperties &rod_properties, Eigen::MatrixXd &forces) {
 int num_vertices = rod.V.rows();
 forces.setZero(num_vertices, 3);
 
 Eigen::VectorXd theta;
 compute_material_frame_angles(rod, theta);
 
 for (int i = 1; i < num_vertices - 1; ++i) {
 Eigen::Vector3d edge1 = rod.V.row(rod.E(i - 1, 1)) - rod.V.row(rod.E(i - 1, 0));
 Eigen::Vector3d edge2 = rod.V.row(rod.E(i, 1)) - rod.V.row(rod.E(i, 0));
 
 // Compute bending force
 double bending_angle = std::acos(edge1.normalized().dot(edge2.normalized()));
 if (bending_angle > rod_properties.max_bending_angle) {
 Eigen::Vector3d bending_force = rod_properties.bending_stiffness * (bending_angle - rod_properties.max_bending_angle) * edge1.cross(edge2).normalized();
 forces.row(i - 1) += bending_force;
 forces.row(i + 1) += bending_force;
 forces.row(i) -= 2 * bending_force;
 }
 
 // Compute twisting force
 if (i > 1 && i < num_vertices - 2) {
 double twisting_angle = theta(i - 1);
 if (std::abs(twisting_angle) > rod_properties.max_twisting_angle) {
 Eigen::Vector3d twisting_force = rod_properties.twisting_stiffness * (twisting_angle - rod_properties.max_twisting_angle) * edge1.normalized();
 forces.row(i - 1) += twisting_force;
 forces.row(i + 1) += twisting_force;
 forces.row(i) -= 2 * twisting_force;
 }
 }
 }
 }
 */
/**
 void update_rod_configuration(Rod &rod, const RodProperties &rod_properties, const Eigen::MatrixXd &forces) {
 int num_vertices = rod.V.rows();
 for (int i = 0; i < num_vertices-1; i++) {
 Eigen::RowVector3d delta = rod_properties.timestep * forces.row(i);
 if(i==0){
 rod.V.row(i) += delta;
 continue;
 }
 rod.V.row(i) += delta;
 // Enforce constraint to maintain edge lengths
 Eigen::RowVector3d prev_edge = rod.V.row(i) - rod.V.row(i - 1);
 Eigen::RowVector3d next_edge = rod.V.row(i + 1) - rod.V.row(i);
 
 prev_edge.normalize();
 next_edge.normalize();
 
 double prev_length = rod.L(i - 1);
 double next_length = rod.L(i);
 
 rod.V.row(i) = rod.V.row(i - 1) + prev_length * prev_edge;
 rod.V.row(i + 1) = rod.V.row(i) + next_length * next_edge;
 }
 }
 */
void draw_arrow(const Eigen::RowVector3d &start, const Eigen::RowVector3d &direction, const Eigen::RowVector3d &color, igl::opengl::glfw::Viewer &viewer, double scale = 0.2) {
    Eigen::RowVector3d end = start + scale * direction;
    viewer.data().add_edges(start, end, color);
}

int main(int argc, char *argv[]) {
    Rod rod;
    RodProperties rod_properties = {
        0.1, // bending_stiffness
        0.1, // twisting_stiffness
    };
    // Read vertex positions from the text file
    //To use the program, you need to modify to your absolute path
    read_vertex_positions("/Users/zhaoyanpeng/582final/vertices.txt", rod.V);
    // Create edges between the consecutive vertices
    int num_vertices = rod.V.rows();
    rod.E.resize(num_vertices - 1, 2);
    rod.L.resize(num_vertices - 1);
    for (int i = 0; i < num_vertices - 1; ++i) {
        rod.E(i, 0) = i;
        rod.E(i, 1) = i + 1;
        rod.L(i) = (rod.V.row(i + 1) - rod.V.row(i)).norm();
    }
    
    // Initialize the viewer
    igl::opengl::glfw::Viewer viewer;
    std::vector<Eigen::Matrix3d> bishop_frames;
    Eigen::VectorXd theta;
    theta.setZero(num_vertices-1);
    theta(0) = M_PI / 2;
    theta(num_vertices-2) = -M_PI / 2;
    std::vector<Eigen::Matrix3d> material_frames;
    Eigen::MatrixXd bend;
    std::vector<double> twist;
    compute_forces(rod, rod_properties, bend, twist);
    
    compute_bishop_frames(rod, bishop_frames);
    
    // Visualize the updated rod configuration
    
    viewer.data().add_points(rod.V, Eigen::RowVector3d(0, 0, 0));
    for (unsigned i = 0; i < rod.E.rows(); ++i) {
        viewer.data().add_edges(
                                rod.V.row(rod.E(i, 0)),
                                rod.V.row(rod.E(i, 1)),
                                Eigen::RowVector3d(1, 0, 0)
                                );
    }
    for(unsigned i = 0; i < rod.E.rows(); ++i){
        Eigen::RowVector3d last_edge_start = rod.V.row(i);
        draw_arrow(last_edge_start, bishop_frames[i].col(0), Eigen::RowVector3d(1, 0, 0), viewer);
        draw_arrow(last_edge_start, bishop_frames[i].col(1), Eigen::RowVector3d(0, 1, 0), viewer);
        draw_arrow(last_edge_start, bishop_frames[i].col(2), Eigen::RowVector3d(0, 0, 1), viewer);
    }
    for(unsigned i = 1; i < rod.E.rows(); ++i){
        Eigen::RowVector3d last_edge_start = rod.V.row(i);
        draw_arrow(last_edge_start, bend.row(i), Eigen::RowVector3d(0, 0, 0), viewer,10);
    }
    for (const auto &element : twist) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
    viewer.launch();
    return 0;
}

