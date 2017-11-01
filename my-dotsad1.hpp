//#ifndef STAN_UDF_FWD_MAT_FUN_DOTPROD_USER_HPP
//#define STAN_UDF_FWD_MAT_FUN_DOTPROD_USER_HPP
/*
#include <stan/math/rev/mat.hpp> // apparently included with stan/math.hpp
#include <Eigen/Dense>
#include <vector>

using Eigen::Matrix;
using Eigen::Dynamic;
*/

    inline double mydotprod(double x, double y, std::ostream* pstream) {
        double out;
        out = x * y;
        return out;
    }

    inline var mydotprod(const var& v1, const var& v2, std::ostream* pstream) {
        // Convert row_vector, vector -> Eigen matrices
        double   x;
        double   y;
        x = v1.val();
        y = v2.val();

        double m = mydotprod(x, y, pstream);
        double dm_dx = v2.val();
        double dm_dy = v1.val();

        return var(new precomp_vv_vari(
                m,                                // value of the output
                v1.vi_,                           // input gradient wrt v1
                v2.vi_,                           // input gradient wrt v2
                dm_dx,                            // partial introduced by this fn wrt v1,
                dm_dy                             // partial introduced by this fn wrt v2
            ));
    }

//#endif