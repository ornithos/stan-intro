import numpy as np

def bspline_deboor(x, t, c, k):
    # ----------------------------
    # k is order not degree!!!!!
    # ----------------------------
    nknots = len(t)
    assert all(t[i] <= t[i+1] for i in range(nknots-1))
    b = np.zeros((k, 1))
    b[0] = 1

    # find (first) knot that is less than x
    left = 0
    assert x > t[0], "query point is not within knot range"
    assert x < t[-1], "query point is not within knot range"

    while x > t[left+1]:
        left += 1

    deltaR = np.zeros((k,1))
    deltaL = np.zeros((k,1))

    for r in range(k-1):
        j         = r + 1
        deltaR[j-1] = t[left+1 + j -1] - x
        deltaL[j-1] = x - t[left+1 + 1 - j -1]
        saved = 0
        for s in range(r+1):
            i     = s + 1
            term  = b[i -1]/(deltaR[i -1] + deltaL[j+1-i -1])
            b[i -1]  = saved + deltaR[i -1]*term
            saved = deltaL[j+1-i -1]*term
        b[j+1 -1] = saved

    out = np.zeros(len(t)-(k-1)+1, dtype=b.dtype)
    out[(left-(k-1)+1):left+2] = b.ravel()
    fnval = sum( (b[i] * c[left-k+2+i] for i in range(k)) )
    return out, fnval