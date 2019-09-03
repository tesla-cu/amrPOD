import numpy as np

def compute_H(X, Psi, i, idx, n, nt, clvl, finest, d_l, G_mat, nl, H):
    
    if clvl < finest-1:
        for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):

            if nl[clvl, idx] > 0:

                l_sum = 0
                for m in range(nl[clvl, idx]):
                    k = G_mat[clvl, idx, m]
                    l_sum += X[j,k] * Psi[k,n]

                l_idx = int([j-i]/d_l[clvl])
                for m in range(l_idx*d_l[clvl], (l_idx+1)*d_l[clvl]):
                    H[m, clvl] = l_sum

            nccells = 0
            for l in range(clvl+1):
                nccells += nl[l, idx]

            if nccells < nt:
                H = compute_H(X, Psi, j, idx, n, nt, clvl+1, finest, d_l, G_mat, nl, H)
            idx += d_l[clvl + 1]

# ------- Good until here
            
    else:

        for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):

            if nl[clvl, idx] > 0:
                l_sum = 0
                for m in range(nl[clvl, idx]):

                    k = G_mat[clvl, idx, m]
                    l_sum += X[j,k] * Psi[k,n]

                for m in range((idx)*d_l[clvl], (idx+1)*d_l[clvl]):
                    H[m, clvl] = l_sum

            if nl[finest, idx] > 0:

                for k in range(j, j+d_l[clvl]):
                    l_sum = 0
                    for m in range(nl[finest, idx]):
                        p = G_mat[finest, idx, m]
                        l_sum += X[k,p] * Psi[p,n]

                    H[k-j+d_l[finest-1] * (idx), finest] = l_sum

            idx += 1

    return H

