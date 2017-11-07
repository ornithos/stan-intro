#include <vector>


void spline_eval(double *t, double x, int T, int k, int m, int *left, double *basis) {
    /* --------------------------------------------
     * t      - (array) knots
     * x      - (double) evaluation point
     * T      - (int) num of knots
     * k      - (int) degree (n.b. not order!)
     * m      - (int) derivative
     * left   - (int ptr) first non-zero basis vector
     * basis  - (array) result / output
     * --------------------------------------------
     */
    double hh[k+1];
    double *h = basis;
    double xb, xa, w;
    int ind, j, n;

    /* find corresponding knot region t[left] <= x <= t[left+1]*/
    *left = 0;
    for (int ii = 0; ii < T-1; ++ii) {
        if (x > t[ii+1]) {
            *left += 1;
            continue;
        }
        /* a little obscure: if equality, keep going, but if
         * a larger value is found, place ptr just before it.
         * This ensures a non-empty region if possible. */
        if (x < t[ii+1]) {
            if (x == t[ii]) {
                *left = ii;
            }
            break;
        }
    }

    /*
     * From scipy. A little different to deBoor (1978/2001):
     * deBoor uses only a single array h, *but* uses a vector for
     * xa and xb, so kind of evens out in favour of this version. (I think.)
     * -------------------------------------- scipy:
     * Perform k-m "standard" deBoor iterations
     * so that h contains the k+1 non-zero values of beta_{ell,k-m}(x)
     * needed to calculate the remaining derivatives.
     */
    h[0] = 1.0;
    for (j = 1; j <= k - m; j++) {
        memcpy(hh, h, j*sizeof(double));
        h[0] = 0.0;
        for (n = 1; n <= j; n++) {
            ind = *left + n;
            xb = t[ind];
            xa = t[ind - j];
            if (xb == xa) {
                h[n] = 0.0;
                continue;
            }
            w = hh[n - 1]/(xb - xa);
            h[n - 1] += w*(xb - x);
            h[n] = w*(x - xa);
        }
    }

    /*
     * Now do m "derivative" recursions
     * to convert the values of beta into the mth derivative
     */
    for (j = k - m + 1; j <= k; j++) {
        memcpy(hh, h, j*sizeof(double));
        h[0] = 0.0;
        for (n = 1; n <= j; n++) {
            ind = *left + n;
            xb = t[ind];
            xa = t[ind - j];
            if (xb == xa) {
                h[m] = 0.0;
                continue;
            }
            w = j*hh[n - 1]/(xb - xa);
            h[n - 1] -= w;
            h[n] = w;
        }
    }
}

void function_eval(double *x, double *t, int k, double *result,
                    double *coeffs, int T, int N) {

    /* --------------------------------------------
     * x      - (dbl array) evaluation array
     * t      - (dbl array) knots
     * k      - (int) degree (n.b. not order!)
     * result - (dbl array) result array (output)
     * coeffs - (dbl array) coefficient array
     * T      - (int) dim. of t
     * N      - (int) dim. of x (these to be replaced when using vector/eigen type)
     * --------------------------------------------
     */
    double basis[k+1];
    int left      = 0;
    int deriv     = 0;

    /* Actual spline evaluation */
    for(int el = 0; el < N; el++) {

        /* deal with boundaries (constant strategy) */
        double xx = x[el];
        if(x[el] < t[k]) {
            xx = t[k];
        } else if(x[el] > t[T-k]) {
            xx = t[T-k];
        }

        spline_eval(t, xx, T, k, deriv, &left, basis);

        /* calculate fnval  */
        double fnval = 0;
        for(int i=0; i < k + 1; i++) {
            fnval += basis[i] * coeffs[left -k + i];
        }
        result[el] = fnval;
    }
}

/* Evaluate derivative: two differences:
 * (i)  deriv = 1 (obviously)
 * (ii) out-of-bounds strategy = 0, rather than const extrapolation */
void deriv_eval(double *x, double *t, int k, double *result,
                    double *coeffs, int T, int N) {
    /* --------------------------------------------
     * x      - (dbl array) evaluation array
     * t      - (dbl array) knots
     * k      - (int) degree (n.b. not order!)
     * result - (dbl array) result array (output)
     * coeffs - (dbl array) coefficient array
     * T      - (int) dim. of t
     * N      - (int) dim. of x (these to be replaced when using vector/eigen type)
     * --------------------------------------------
     */
    double basis[k+1];
    int left      = 0;
    int deriv     = 1;

    /* Actual spline evaluation */
    for(int el = 0; el < N; el++) {

        /* deal with boundaries (constant strategy)
         * (NOTICE DIFFERENCE TO function_eval) */
        double xx = x[el];
        if(x[el] < t[k]) {
            result[el] = 0;
            continue;
        } else if(x[el] > t[T-k]) {
            result[el] = 0;
            continue;
        }

        spline_eval(t, xx, T, k, deriv, &left, basis);

        /* calculate fnval  */
        double fnval = 0;
        for(int i=0; i < k + 1; i++) {
            fnval += basis[i] * coeffs[left -k + i];
        }
        result[el] = fnval;
    }
}

/* ---- OUT OF BOUNDS NOT CHECKED FOR. COEFF VECTOR MUST BE OF LENGTH T-k, etc. ---- */