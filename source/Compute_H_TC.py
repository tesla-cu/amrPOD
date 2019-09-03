import numpy as np

def compute_H_TC(X, i, idx, n, nt, clvl, finest, d_l, G_mat, nl, H, Phi_count_imp_arith, Phi_count_imp_access, Phi_count_imp_assign, Phi_count_imp_logtest, Phi_count_imp_func):
    
    if clvl < finest-1:
        Phi_count_imp_arith   += 1
        Phi_count_imp_logtest += 1

        for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
            Phi_count_imp_arith  += 1
            Phi_count_imp_access += 1

            if nl[clvl, idx] > 0:
                Phi_count_imp_logtest += 1
                Phi_count_imp_access  += 1

                l_sum = 0
                Phi_count_imp_assign += 1

                for m in range(nl[clvl, idx]):
                    Phi_count_imp_arith  += 1
                    Phi_count_imp_access += 1

                    k = G_mat[clvl, idx, m]
                    Phi_count_imp_access += 1
                    Phi_count_imp_assign += 1

                    #l_sum += X[j,k] * Psi[k,n]
                    Phi_count_imp_access += 2
                    Phi_count_imp_arith  += 2
                    Phi_count_imp_assign += 1

                l_idx = int([j-i]/d_l[clvl])
                Phi_count_imp_arith  += 2
                Phi_count_imp_access += 2
                Phi_count_imp_assign += 1

                for m in range(l_idx*d_l[clvl], (l_idx+1)*d_l[clvl]):
                    Phi_count_imp_arith  += 1
                    Phi_count_imp_access += 1

                    #H[m, clvl] = l_sum
                    Phi_count_imp_access += 1
                    Phi_count_imp_assign += 1

            nccells = 0
            Phi_count_imp_assign += 1

            for l in range(clvl+1):
                Phi_count_imp_arith  += 1
                Phi_count_imp_access += 1

                nccells += nl[l, idx]
                Phi_count_imp_access += 1
                Phi_count_imp_arith  += 1
                Phi_count_imp_assign += 1

            if nccells < nt:
                Phi_count_imp_logtest += 1

                H = compute_H(X, j, idx, n, nt, clvl+1, finest, d_l, G_mat, nl, H, Phi_count_imp_arith, Phi_count_imp_access, Phi_count_imp_assign, Phi_count_imp_logtest, Phi_count_imp_func)
                Phi_count_imp_func += 1

            idx += d_l[clvl + 1]
            Phi_count_imp_arith  += 2
            Phi_count_imp_access += 1
            Phi_count_imp_assign += 1

# ------- Good until here
            
    else:

        for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
            Phi_count_imp_arith  += 1
            Phi_count_imp_access += 1

            if nl[clvl, idx] > 0:
                Phi_count_imp_access  += 1
                Phi_count_imp_logtest += 1

                #l_sum = 0
                Phi_count_imp_assign += 1

                for m in range(nl[clvl, idx]):
                    Phi_count_imp_arith  += 1
                    Phi_count_imp_access += 1

                    k = G_mat[clvl, idx, m]
                    Phi_count_imp_access += 1
                    Phi_count_imp_assign += 1

                    #l_sum += X[j,k] * Psi[k,n]
                    Phi_count_imp_arith  += 2
                    Phi_count_imp_access += 2
                    Phi_count_imp_assign += 1

                for m in range((idx)*d_l[clvl], (idx+1)*d_l[clvl]):
                    Phi_count_imp_arith  += 1
                    Phi_count_imp_access += 1

                    # H[m, clvl] = l_sum
                    Phi_count_imp_access += 1
                    Phi_count_imp_assign += 1

            if nl[finest, idx] > 0:
                Phi_count_imp_access  += 1
                Phi_count_imp_logtest += 1

                for k in range(j, j+d_l[clvl]):
                    Phi_count_imp_arith  += 1
                    Phi_count_imp_access += 1

                    #l_sum = 0
                    Phi_count_imp_assign += 1

                    for m in range(nl[finest, idx]):
                        Phi_count_imp_arith  += 1
                        Phi_count_imp_access += 1

                        p = G_mat[finest, idx, m]
                        Phi_count_imp_access += 1
                        Phi_count_imp_assign += 1

                        #l_sum += X[k,p] * Psi[p,n]
                        Phi_count_imp_access += 2
                        Phi_count_imp_arith  += 2
                        Phi_count_imp_assign += 1

                    # H[k-j+d_l[finest-1] * (idx), finest] = l_sum
                    Phi_count_imp_access += 2
                    Phi_count_imp_arith  += 4
                    Phi_count_imp_assign += 1

            idx += 1
            Phi_count_imp_arith += 1

    return H, Phi_count_imp_arith, Phi_count_imp_access, Phi_count_imp_assign, Phi_count_imp_logtest, Phi_count_imp_func

