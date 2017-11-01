#ifndef STAN_UDF_FWD_MAT_FUN_DOTPROD_USER_HPP
#define STAN_UDF_FWD_MAT_FUN_DOTPROD_USER_HPP

#include <stan/math/rev/mat.hpp> // apparently included with stan/math.hpp
#include <Eigen/Dense>
#include <vector>

using Eigen::Matrix;
using Eigen::Dynamic;

namespace stan {
  namespace math {

    inline double mydotprod(const Matrix<double, 1, Dynamic>& x, const Matrix<double, Dynamic, 1>& y, std::ostream* pstream) {
        double out;
        out = x * y;
        return out;
    }

    inline var mydotprod(const Matrix<var, 1, Dynamic>& v1, const Matrix<var, Dynamic, 1>& v2, std::ostream* pstream) {
        // Convert row_vector, vector -> Eigen matrices
        Matrix<double, 1, Dynamic>   x;
        Matrix<double, Dynamic, 1>   y;
        x = value_of(v1);
        y = value_of(v2);
    //    Matrix<double, 1, Dynamic>   x(v1.size());
    //    Matrix<double, Dynamic, 1>   y(v2.size());
    //    for (int i=0; i < x.size(); ++i) x[i] = v1[i].val();
    //    for (int i=0; i < x.size(); ++i) y[i] = v2[i].val();
    //        Matrix<double, 1, Dynamic> x;
    //        Matrix<double, Dynamic, 1> y;
    //        x = value_of(v1);
    //        y = value_of(v2);
        // call function
        double m = mydotprod(x, y, pstream);

        // partial derivatives
    //    Matrix<double, 1, Dynamic>        dm_dx = y;
    //    Matrix<double, 1, Dynamic>        dm_dy = x.transpose();

        // autodiff wrapper
    //    std::vector<Matrix<double, 1, Dynamic> > grads(2);
    //    grads[0] = dm_dx;
    //    grads[1] = dm_dy;
    //    std::vector<std::vector<var> > v(2);
    //    v[0] = v1;
    //    v[1] = v2;


    //    std::vector<var>
    //    MatrixXd V(2*y.rows(),1);
    //    V << x.transpose(), y;
    //    MatrixXd D(2*dm_dy.cols(),1);
    //    D << dm_dx.transpose(), dm_dy.transpose();

        // need to copy as v1 is const

    //    std::vector<Matrix<var, Dynamic, 1> > V(2);
    //    V[0] = v1;
    //    V[1] = v2;
        std::vector<var> V;
        for(int i = 0; i < x.rows(); ++i)
            V.push_back(v1(1,i));
        for(int i = 0; i < x.rows(); ++i)  // fine to be x since x,y same len
            V.push_back(v2(i,1));

        std::vector<double> D;
        for(int i = 0; i < x.rows(); ++i)
            D.push_back(y(i,1));
        for(int i = 0; i < x.rows(); ++i)
            D.push_back(x(1,i));
    //
    //
    //    std::vector<double> Dx(y.data(), y.data() + y.rows() * y.cols());
    //    std::vector<double> Dy(x.data(), x.data() + x.rows() * x.cols());
    //    Dx.insert(Dx.end(), Dy.begin(), Dy.end());


        var retval;
    //    retval = precomputed_gradients(m, V, D);
        retval = precomputed_gradients(m, V, D);
        return retval;
    //    return var(new precomp_vv_vari(
    //        m,                                // value of the output
    //        x,                                // input gradient wrt v1
    //        y,                                // input gradient wrt v2
    //        dm_dx,                            // partial introduced by this fn wrt v1,
    //        dm_dy                             // partial introduced by this fn wrt v2
    //    ));
    }
  }
}
#endif