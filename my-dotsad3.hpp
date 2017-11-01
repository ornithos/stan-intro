//#ifndef STAN_UDF_FWD_MAT_FUN_DOTPROD_USER_HPP
//#define STAN_UDF_FWD_MAT_FUN_DOTPROD_USER_HPP

#include <stan/math/rev/mat.hpp> // apparently included with stan/math.hpp
#include <Eigen/Dense>
#include <vector>

using Eigen::Matrix;
using Eigen::Dynamic;


    inline double mydotprod(const Matrix<double, 1, 1>& x, const Matrix<double, 1, 1>& y, std::ostream* pstream) {
        std::cout << "using double double version" << std::endl;
        double out;
        out = x * y;
        return out;
    }

    inline var mydotprod(const Matrix<var, 1, 1>& v1, const Matrix<var, 1, 1>& v2, std::ostream* pstream) {
        // Convert row_vector, vector -> Eigen matrices
        std::cout << "using var var version" << std::endl;
        Matrix<double, 1, 1>   x;
        Matrix<double, 1, 1>   y;
        x = value_of(v1);
        y = value_of(v2);

        double m = mydotprod(x, y, pstream);
        double dm_dx = y(0,0);
        double dm_dy = x(0,0);

        return var(new precomp_vv_vari(
                m,                                // value of the output
                v1(0,0).vi_,                           // input gradient wrt v1
                v2(0,0).vi_,                           // input gradient wrt v2
                dm_dx,                            // partial introduced by this fn wrt v1,
                dm_dy                             // partial introduced by this fn wrt v2
            ));
    }

    inline var mydotprod(const Matrix<var, 1, 1>& v1, const Matrix<double, 1, 1>& v2, std::ostream* pstream) {
        // Convert row_vector, vector -> Eigen matrices
        std::cout << "using var double version" << std::endl;
        Matrix<double, 1, 1>   x;
        x = value_of(v1);

        double m = mydotprod(x, v2, pstream);
        double dm_dx = v2(0,0);

        return var(new precomp_v_vari(
                m,                                // value of the output
                v1(0,0).vi_,                           // input gradient wrt v1
                dm_dx                            // partial introduced by this fn wrt v1,
            ));
    }

    inline var mydotprod(const Matrix<double, 1, 1>& v1, const Matrix<var, 1, 1>& v2, std::ostream* pstream) {
        // Convert row_vector, vector -> Eigen matrices
        std::cout << "using double var version" << std::endl;
        Matrix<double, 1, 1>   y;
        y = value_of(v2);

        double m = mydotprod(v1, y, pstream);
        double dm_dy = v1(0,0);

        return var(new precomp_v_vari(
                m,                                // value of the output
                v2(0,0).vi_,                           // input gradient wrt v2
                dm_dy                            // partial introduced by this fn wrt v2,
            ));
    }

//#endif