#include <iostream>
#include <vector>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/gaussian_curvature.h>
#include <igl/opengl/glfw/Viewer.h>

#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

using namespace std;
using namespace Eigen;
using namespace igl;

void compute_gravity_center(Vector3d& p0, Vector3d& p1, Vector3d& p2,
                            Vector3d& gravity_center, Vector3d& area_ratio, double& area){
    double l0 = (p1 - p2).norm();//p1 - p2
    double l1 = (p0 - p2).norm();//p0 - p2
    double l2 = (p0 - p1).norm();//p0 - p1

    double mid = (l0 + l1 + l2) / 2;
    area = mid * (mid - l0) * (mid - l1) * (mid - l2);

    double r0 = (pow(l0, 4) - pow(l1 * l1 - l2 * l2, 2) ) / (16 * area);
    double r1 = (pow(l1, 4) - pow(l2 * l2 - l0 * l0, 2)) / (16 * area);
    double r2 = (pow(l2, 4) - pow(l0 * l0 - l1 * l1, 2)) / (16 * area);
    gravity_center = r0 * p0 + r1 * p1 + r2 * p2;
    area_ratio << r0, r1, r2;
}

int main() {

    MatrixXd V;
    MatrixXi F;
    readOFF("../data/bumpy.off", V, F);
    MatrixXd local_area_v(F.rows(), 3); // #face_num, 3
    for(int i = 0 ; i < F.rows() ; i ++){
        //ith facet
        Vector3d p0 = V.row(F(i, 0));
        Vector3d p1 = V.row(F(i, 1));
        Vector3d p2= V.row(F(i, 2));
        double area;
        Vector3d gravity_center, area_ratio;
        compute_gravity_center(p0, p1, p2, gravity_center, area_ratio, area);
        area = sqrt(area);
        bool flag = true;
        for(int k = 0 ; k < area_ratio.size() ; k ++){
            if(area_ratio(k) < 0) flag = false;
        }
        for(int j = 0 ; j < 3 ;j ++){
            if(flag)
                local_area_v(i, j) = area * (1 - area_ratio(j)) / 2;
            else
                local_area_v(i,j) = area * 0.5;
        }
    }
//    cout << local_area_v<< endl;
    vector<vector<int>> VF;
    vector<vector<int>> VFI;
    vertex_triangle_adjacency( V, F, VF, VFI);
    VectorXd gassian_curvatures(V.rows(), 1);
    for(int i = 0 ; i < V.rows() ; i ++){
        Vector3d p0 = V.row(i);
        vector<int> adjacency_facets = VF[i];
        vector<int> adjacency_idxes = VFI[i];

        double angle_sum = 0.0;
        double area_sum = 0.0;
        for(int idx = 0; idx < adjacency_facets.size() ; idx ++){
            Vector3i facet = F.row(adjacency_facets[idx]);
            int vertex_idx_facet = adjacency_idxes[idx];
            area_sum += local_area_v(adjacency_facets[idx], vertex_idx_facet);
            int one, two;
            if(vertex_idx_facet == 0){
                one = 1; two = 2;
            }else if(vertex_idx_facet == 1){
                one = 0; two = 2;
            }else{
                one = 0; two = 1;
            }
            Vector3d p1 = V.row(facet(one));
            Vector3d p2 = V.row(facet(two));
            Vector3d p01 = p1 - p0;
            Vector3d p02 = p2 - p0;

            double theta = acos ((p01.dot(p02) ) / (p01.norm() * p02.norm()));
            angle_sum += theta;


        }
        gassian_curvatures(i) = (2 * PI - angle_sum) ;
    }
//    cout << gassian_curvatures << endl;
    gassian_curvatures = gassian_curvatures.normalized();
    VectorXd K;
//    gaussian_curvature(V, F, K);
    K = gassian_curvatures;
    SparseMatrix<double>  M,Minv;
    massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
    invert_diag(M,Minv);
    // Divide by area to get integral average
    K = (Minv*K).eval();

    opengl::glfw::Viewer viewer;

    viewer.data().set_mesh(V, F);
    viewer.data().set_data(K);
    viewer.launch();
//    cout << K <<endl;
    return 0;
}
