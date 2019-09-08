import numpy as np

def compute_A_TC(X_grid, R, d_l, nt, nspat, finest):

        A_count_unalt        = 0
        A_count_unalt_arith  = 0
        A_count_unalt_access = 0
        A_count_unalt_assign = 0

        A_count_imp          = 0
        A_count_imp_arith    = 0
        A_count_imp_access   = 0
        A_count_imp_assign   = 0
        A_count_imp_logtest  = 0        

        #--------- Unaltered 


        #A_unalt = np.zeros((nt, nt))
        A_count_unalt_assign += 1

        for i in range(nt):
                A_count_unalt_arith  += 1
                A_count_unalt_access += 1

                for j in range(nt):
                        A_count_unalt_arith  += 1
                        A_count_unalt_access += 1

                        #a_sum = 0
                        A_count_unalt_assign += 1

                        for k in range(nspat):
                                A_count_unalt_arith  += 1
                                A_count_unalt_access += 1

                                #a_sum += X[k,i] * Phi[k,j]
                                A_count_unalt_arith  += 2
                                A_count_unalt_access += 2
                                A_count_unalt_assign += 1

                        #A_unalt[i,j] = a_sum
                        A_count_unalt_access += 1
                        A_count_unalt_assign += 1


#--------- Implemented
        jloop_count     = 0
        if_X_grid_count = 0
        Precomp_count   = 0
        #A_imp = np.zeros((nt,nt))
        A_count_imp_assign += 1
        G = np.zeros((nspat), dtype=int)
        A_count_imp_assign += 1

        i     = 0
        A_count_imp_assign += 1

        for n in range(nspat):
                A_count_imp_arith  += 1
                A_count_imp_access += 1

                A_count_imp_logtest +=1
                if i < nspat:

                        X_grid_max = X_grid[i,0]
                        A_count_imp_access += 1
                        A_count_imp_assign += 1

                        A_count_imp_logtest +=1
                        if X_grid_max < finest:

                                for j in range(1,nt):
                                        jloop_count += 1
                                        A_count_imp_arith  += 1
                                        A_count_imp_access += 1

                                        A_count_imp_logtest  += 1
                                        A_count_imp_access   += 1       
                                        if X_grid[i,j] > X_grid_max:
                                                if_X_grid_count += 1

                                                X_grid_max = X_grid[i,j]
                                                A_count_imp_access += 1
                                                A_count_imp_assign += 1

                                                A_count_imp_logtest +=1
                                                if X_grid_max == finest:

                                                        break

                        G[i]  = d_l[X_grid_max]
                        A_count_imp_access += 2
                        A_count_imp_assign += 1

                        i    += G[i]
                        A_count_imp_access += 1
                        A_count_imp_assign += 1
                        A_count_imp_arith  += 1
                else:
                        break

        Precomp_count += A_count_imp_access + A_count_imp_arith + A_count_imp_assign + A_count_imp_logtest
        # print('Precomp ops = ', Precomp_count)



        # for n in range(nspat):
        #       A_count_imp_arith  += 1
        #       A_count_imp_access += 1

        #       A_count_imp_logtest +=1
        #       if i < nspat:

        #               X_grid_max = X_grid[i,0]
        #               A_count_imp_access += 1
        #               A_count_imp_assign += 1

        #               for j in range(1,nt):
        #                       A_count_imp_arith  += 1
        #                       A_count_imp_access += 1

        #                       A_count_imp_logtest  += 1
        #                       A_count_imp_access += 1 
        #                       if X_grid[i,j] > X_grid_max:                            

        #                               X_grid_max = X_grid[i,j]
        #                               A_count_imp_access += 1
        #                               A_count_imp_assign += 1

        #               G[i]  = d_l[X_grid_max]
        #               A_count_imp_access += 2
        #               A_count_imp_assign += 1

        #               i    += G[i]
        #               A_count_imp_access += 1
        #               A_count_imp_assign += 1
        #               A_count_imp_arith  += 1

        #       else:
        #               break

        for i in range(nt):
                A_count_imp_arith  += 1
                A_count_imp_access += 1

                for j in range(nt):
                        A_count_imp_arith  += 1
                        A_count_imp_access += 1

                        #a_sum = 0
                        A_count_imp_assign += 1

                        k     = 0
                        A_count_imp_assign += 1

                        for n in range(nspat):
                                A_count_imp_arith  += 1
                                A_count_imp_access += 1

                                A_count_imp_logtest += 1
                                if k < nspat:

                                        #a_sum += G[k] * X[k,i] * Phi[k,j]
                                        A_count_imp_arith  += 3
                                        A_count_imp_access += 3
                                        A_count_imp_assign += 1

                                        k += G[k]
                                        A_count_imp_access += 1
                                        A_count_imp_arith  += 1
                                        A_count_imp_assign += 1



                                else:
                                        break
                        #A_imp[i,j] = a_sum
                        A_count_imp_access += 1
                        A_count_imp_assign += 1

        A_count_unalt = A_count_unalt_access + A_count_unalt_assign + A_count_unalt_arith
        A_count_imp   = A_count_imp_access + A_count_imp_assign + A_count_imp_arith + A_count_imp_logtest

        # print('jloop_count = ', jloop_count)
        # print('if_X_grid_count = ', if_X_grid_count)

        return A_count_imp, A_count_unalt 
