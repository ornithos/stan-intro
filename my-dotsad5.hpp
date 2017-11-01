//#ifndef STAN_UDF_FWD_MAT_FUN_DOTPROD_USER_HPP
//#define STAN_UDF_FWD_MAT_FUN_DOTPROD_USER_HPP

//#define ALEX_UDF_VERBOSE_LEVEL = True;

/*
 * THIS APPEARS TO BE SUBSTANTIALLY SLOWER THAN USING THE PRECOMP_VV_VARI FORMULATION
 * ----------> !!!!!!!!!!!!!! Don't use !!!!!!!!!!!!!!!!!!!!!!!!! <------------------
 */

 /*
  * NOT SURE WHAT TO MAKE OF THE SPEED (sampling), BUT BOTH BUG OUT WHEN USING OPTIMIZE.
  * -----------------------------------------------------------------------------------
  * My current observation is that the various signatures still use the same return type
  * in precomp grads. Worth looking at the `decouple_ode_states` example -- uses (perhaps?)
  * old version of boost:promote_args for returning the vars (i.e. grad is always type
  * double, but node type changes.)
  */

#include <stan/math/rev/mat.hpp> // apparently included with stan/math.hpp
#include <Eigen/Dense>
#include <vector>

using Eigen::Matrix;
using Eigen::Dynamic;
using boost::math::tools::promote_args;


    inline double mydotprod(const Matrix<double, 1, Dynamic>& x, const Matrix<double, Dynamic, 1>& y, std::ostream* pstream) {
        //std::cout << "using double double version" << std::endl;
        double out;
        out = x * y;
        return out;
    }

    inline var mydotprod(const Matrix<var, 1, Dynamic>& v1, const Matrix<var, Dynamic, 1>& v2, std::ostream* pstream) {
        // Convert row_vector, vector -> Eigen matrices
        //std::cout << "using var var version" << std::endl;
        Matrix<double, 1, Dynamic>   x;
        Matrix<double, Dynamic, 1>   y;
        x = value_of(v1);
        y = value_of(v2);

        double m = mydotprod(x, y, pstream);
        double l = x.cols();

        std::vector<var> V;
        for(int i = 0; i < l; ++i)
            V.push_back(v1(1,i));
        for(int i = 0; i < l; ++i)
            V.push_back(v2(i,1));

        std::vector<double> D;
        for(int i = 0; i < l; ++i)
            D.push_back(y(i,1));
        for(int i = 0; i < l; ++i)
            D.push_back(x(1,i));

        var retval;
        retval = precomputed_gradients(m, V, D);
        return retval;
    }

    inline var mydotprod(const Matrix<var, 1, Dynamic>& v1, const Matrix<double, Dynamic, 1>& v2, std::ostream* pstream) {
        // Convert row_vector, vector -> Eigen matrices
        std::cout << "using var double version" << std::endl;
        Matrix<double, 1, Dynamic>   x;
        Matrix<double, Dynamic, 1>   y;
        x = value_of(v1);
        y = value_of(v2);

        double m = mydotprod(x, v2, pstream);

        std::vector<var> V;
        for(int i = 0; i < x.rows(); ++i)
            V.push_back(v1(1,i));

        std::vector<double> D;
        for(int i = 0; i < x.rows(); ++i)
            D.push_back(y(1,i));

        var retval;
        retval = precomputed_gradients(m, V, D);
        return retval;
    }

    inline var mydotprod(const Matrix<double, 1, Dynamic>& v1, const Matrix<var, Dynamic, 1>& v2, std::ostream* pstream) {
        // Convert row_vector, vector -> Eigen matrices
        std::cout << "using double var version" << std::endl;
        Matrix<double, 1, Dynamic>   x;
        Matrix<double, Dynamic, 1>   y;
        x = value_of(v1);
        y = value_of(v2);

        double m = mydotprod(v1, y, pstream);

        std::vector<var> V;
        for(int i = 0; i < x.rows(); ++i)
            V.push_back(v2(1,i));

        std::vector<double> D;
        for(int i = 0; i < x.rows(); ++i)
            D.push_back(x(1,i));

        var retval;
        retval = precomputed_gradients(m, V, D);
        return retval;
    }

//#endif