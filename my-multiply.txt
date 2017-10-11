inline double mymultiply(double x, double y, std::ostream* pstream) {
    double out;
    out = x * y;
    return out;
}

inline var mymultiply(const var& v1, const var& v2, std::ostream* pstream) {
    double x = v1.val(),
           y = v2.val(),
           m = mymultiply(x, y, pstream);

    // partial derivatives
    double dm_dx = y,
           dm_dy = x;

    // autodiff wrapper
    return var(new precomp_vv_vari(
        m,                                // value of the output
        v1.vi_,                           // input gradient wrt v1
        v2.vi_,                           // input gradient wrt v2
        dm_dx,                            // partial introduced by this fn wrt v1,
        dm_dy                             // partial introduced by this fn wrt v2
    ));
}