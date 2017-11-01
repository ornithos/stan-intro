#ifndef STAN_UDF_FWD_MAT_FUN_DOTPROD_USER_HPP
#define STAN_UDF_FWD_MAT_FUN_DOTPROD_USER_HPP

#include <stan/math/rev/mat.hpp> // apparently included with stan/math.hpp
#include <Eigen/Dense>
#include <vector>

using Eigen::Matrix;
using Eigen::Dynamic;

namespace stan {
  namespace math {

    inline double mydotprod(const Matrix<double, 1, 1>& x, const Matrix<double, 1, 1>& y, std::ostream* pstream) {
        double out;
        out = x * y;
        return out;
    }

    inline var mydotprod(const Matrix<var, 1, 1>& v1, const Matrix<var, 1, 1>& v2, std::ostream* pstream) {
        // Convert row_vector, vector -> Eigen matrices
        Matrix<double, 1, 1>   x;
        Matrix<double, 1, 1>   y;
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

        return var(new precomp_vv_vari(
                m,                                // value of the output
                v1,                           // input gradient wrt v1
                v2,                           // input gradient wrt v2
                y,                            // partial introduced by this fn wrt v1,
                x                             // partial introduced by this fn wrt v2
            ));
    }
  }
}
#endif