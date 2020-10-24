#include <iostream>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
using namespace std;
using namespace Eigen;
using namespace igl;


Vector3d compute_gravity_center(MatrixXd& points){
    Matrix3d D(3, 3);
    Vector3d p;

    double A1 = 2 * (points(1, 0) - points(0, 0)),
           B1 = 2 * (points(1, 1) - points(0, 1)),
           C1 = 2 * (points(1, 2) - points(0, 2)),
           D1 = points(1, 0) * points(1, 0) + points(1, 1) * points(1, 1) + points(1, 2) * points(1 , 2)
                   - points(0, 0) * points(0, 0) - points(0, 1) * points(0, 1) - points(0, 2) * points(0, 2);
    double A2 = 2 * (points(2, 0) - points(0, 0)),
           B2 = 2 * (points(2, 1) - points(0, 1)),
           C2 = 2 * (points(2, 2) - points(0, 2)),
           D2 = points(2, 0) * points(2, 0) + points(2, 1) * points(2, 1) + points(2, 2) * points(2, 2)
                   - points(0, 0) * points(0, 0) - points(0, 1) * points(0, 1) - points(0, 2) * points(0, 2);
    double A3 = 2 * (points(2, 0) - points(1, 0)),
           B3 = 2 * (points(2, 1) - points(1, 1)),
           C3 = 2 * (points(2, 2) - points(1, 2)),
           D3 = points(2, 0) * points(2, 0) + points(2, 1) * points(2, 1) + points(2, 2) * points(2, 2)
                   - points(1, 0) * points(1, 0) - points(1, 1) * points(1, 1) - points(1, 2) * points(1 , 2);


    Vector3d p01 = points.row(1) - points.row(0);
    Vector3d p02 = points.row(2) - points.row(0);


    D(0, 0) = A1;D(0, 1) = B1;D(0, 2) = C1;
    D(1, 0) = A2; D(1, 1) = B2; D(1, 2) = C2;
    D(2, 0) = A3;D(2, 1) = B3; D(2, 2) = C3;
    cout << "D is " << D <<endl;
    cout << D.determinant() << endl;


    Vector3d right = {D1, D2, D3};
    cout << "right : " << right << endl;
//    Vector3d answer = D.colPivHouseholderQr().solve(right);
    Vector3d answer = D.llt().solve(right);
    cout << "Answer" << answer << endl;
    if(D.determinant() != 0){

        for(int i = 0 ; i < 3 ; i ++){
            cout << "##### i " << i << endl;
            Matrix3d D_i = D;
            D_i(0, i) = D1;D_i(1, i) = D2;D_i(2, i) = D3;
            cout << "Current D_i" << D_i <<endl;
            p[i] = D_i.determinant() / D.determinant();
        }
    }

    return p;

}

int main() {

    MatrixXd V;
    MatrixXi F;
    MatrixXd local_area_v;
    readOBJ("../data/elephant.obj", V, F);



    for(int i = 0 ; i < 1 ; i ++){
        //ith facet
        MatrixXd M(3, 3);
//        M.row(0) = V.row(F(i, 0));
//        M.row(1) = V.row(F(i, 1));
//        M.row(2) = V.row(F(i, 2));
        M.row(0) = Vector3d{1., 0., 0.};
        M.row(1) = Vector3d{0., 1., 0.};
        M.row(2) = Vector3d(0., 0., 1.);
        cout << M << endl;

        Vector3d gcenter = compute_gravity_center(M);
        cout << "垂心" << gcenter << endl;
    }



    cout << endl;
    return 0;
}
