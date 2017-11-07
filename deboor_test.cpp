#define PRINT_TO_CONSOLE 0

#include <stdio.h>      /* printf */
#include <stdlib.h>     /* atoi */
#include <iostream>
#include <vector>
#include <string.h>     /* strcmp*/
#include "deboor_spline.hpp"


int main(int argc, char *argv[]) {

    if(argc == 1) {   // only function name
        printf("Insufficient arguments passed in!\n");
        return -1;
    }

    int k = atoi(argv[1]);
    int T, N;

    // T and k
    if (argc == 2) {
        T = 13;
        N = 6;
    } else {
        if(argc < 4) {   // only function name
            printf("args: k T N [knots] [coeffs] [query]\n");
            return -1;
        }
        T = atoi(argv[2]);
        N = atoi(argv[3]);
    }

    // find if derivative
    const char* char_d = "d";
    int deriv = 0;
    if(strcmp(argv[argc-1], char_d) == 0) {
//        printf("Got derivative flag!\n");
        deriv = 1;
        argc--;
    }

    // define all inputs
    std::vector<double> arrKnots(T);
    std::vector<double> coeffs(T-k);
    double x[N];

    /* initialise vectors: knots and coeffs */
    if (argc == 2) {
        // defaults
        const int T_init              = 13;
        double arrKnots_init[T_init]  = {0,0,0,0,1,2,3,4,5,6,6,6,6};
        double coeffs_init[10]        = {0,0,0,-1,2,0,-1,0,0,0};

        for(int i = 0; i < T; i++) {
            arrKnots[i] = arrKnots_init[i];
        }
        for(int i = 0; i < 10; i++) {
            coeffs[i] = coeffs_init[i];
        }
        if(T-k > 10) {
            for(int i = 10; i < T-k; i++) {
                coeffs[i] = 0;
            }
        }
        //query point
        const int N_init = 6;
        double x_init[N_init] = {1.5,2.0,2.5,3.0,3.5,4.0}; //{1.5,2.0,3.0,3.5,4.0,4.5};
        for(int i = 0; i < N; i++)
            x[i] = x_init[i];
        if(N > N_init) {
            printf("N inconsistent with N_init. See source.\n");
            return -1;
        }

    }
    else {
        int exp_n_args = 2*T - k + N + 3;
        if(argc -4 != exp_n_args -3) {
            printf("Expecting %d args, got %d\n", exp_n_args - 3, argc - 4);
            return - 1;
        }
        #if PRINT_TO_CONSOLE == 1
            printf("Knot values:\n");
        #endif
        for(int ii = 0; ii < T; ii++) {
            arrKnots[ii] = atof(argv[4+ii]);
            #if PRINT_TO_CONSOLE == 1
                printf("%.4f\n", arrKnots[ii]);
            #endif
        }
        #if PRINT_TO_CONSOLE == 1
            printf("Coeff values:\n");
        #endif
        for(int ii = 0; ii < T-k; ii++) {
            coeffs[ii] = atof(argv[4+T+ii]);
            #if PRINT_TO_CONSOLE == 1
                printf("%.4f\n", coeffs[ii]);
            #endif
        }
        #if PRINT_TO_CONSOLE == 1
            printf("Query values:\n");
        #endif
        for(int ii = 0; ii < N; ii++) {
            x[ii] = atof(argv[4+2*T-k+ii]);
            #if PRINT_TO_CONSOLE == 1
                printf("%.4f\n", x[ii]);
            #endif
        }
        #if PRINT_TO_CONSOLE == 1
            printf("\n");
        #endif

    }

    std::vector<double> basis(k+1);
    double result[k+1];
    if (!deriv) {
        bspline::function_eval(x, &arrKnots[0], k, result, &coeffs[0], T, N);
    }
    else {
        bspline::deriv_eval(x, &arrKnots[0], k, result, &coeffs[0], T, N);
    }

    /* Actual spline evaluation */
    for(int el = 0; el < N; el++) {
        /*spline_eval(&arrKnots[0], x[el], k, deriv, &left, &basis[0]);

        // calculate fnval
        double fnval = 0;
        for(int i=0; i < k + 1; i++) {
            fnval += basis[i] * coeffs[left -k + i];
            //printf("%.3f\t", basis[i]);
        }*/

        printf("%.8f\n", result[el]);
    }

    return 0;
}