function sorted2midpoints(sorted)
    n = length(sorted)
    midpoints = zeros(n + 1)
    Δ = (sorted[1] + sorted[n])/(n+1)
    midpoints[1] = sorted[1] - Δ
    midpoints[n+1] = sorted[n] + Δ
    midpoints[2:n] = 0.5(sorted[2:n] + sorted[1:n-1])
    return midpoints
end

function coarseness_g(
        S::Matrix{<:Real}, 
        w::Vector{<:Integer};
        vervose=false::Bool,
        report=""::String
    )
    
    n = length(w)
    
    # Utils ##################################################################################
    sortperm_x = sortperm(S[:, 1])
    sorted_x = S[sortperm_x, 1]
    midpoints_x = sorted2midpoints(sorted_x)

    sortperm_y = sortperm(S[:, 2])
    sorted_y = S[sortperm_y, 2]
    midpoints_y = sorted2midpoints(sorted_y)

    Disc_pqij(p,q,i,j) = abs(sum(w[intersect(sortperm_x[p:p+i-1],sortperm_y[q:q+j-1])]))
    Count_pqij(p,q,i,j) = length(intersect(sortperm_x[p:p+i-1],sortperm_y[q:q+j-1]))

    # Initialization #########################################################################
    inf = 100
    C = 100*ones(Int,n,n,n,n)
    P = zeros(Int,n,n,n,n)
    D = zeros(Int,n,n,n,n)

    # 0 --> not calculated!
    # 1 --> Discepancy
    # 2 --> Vertical
    # 3 --> Horizontal
    argcuts = zeros(Int,n,n,n,n)

    # 0   --> not calculated! or Discepancy
    # 1:n --> midpoint index of v or h cut
    idxcuts = zeros(Int,n,n,n,n)

    for p = 1:n
        for q = 1:n
            for i = 1:n-p+1
                for j = 1:n-q+1
                    P[p,q,i,j] = Count_pqij(p,q,i,j)
                    D[p,q,i,j] = Disc_pqij(p,q,i,j)
                end
            end
            if P[p,q,1,1] != 0
                C[p,q,1,1] = D[p,q,1,1]
            end
        end
    end
    
    
    for i = 1:n
        for j = 1:n
            if i+j > 2
                for p = 1:n-i+1
                    for q = 1:n-j+1
                        # CREO QUE ESTE IF NO CAMBIA NADA; PERO EVITA HARTO CALCULO TRIVIAL
                        if P[p,q,i,j] != 0

                            # vertical cut
                            if i > 1
                                v_cuts = [[C[p,q,s,j], C[p+s,q,i-s,j]][k] for s=1:i-1, k=1:2]
                                for s=1:i-1
                                    if v_cuts[s,1]==inf || v_cuts[s,2]==inf
                                        v_cuts[s,:] = [-1, -1]
                                    end
                                end
                                v_argmins = argmin(v_cuts, dims=2)
                                v_mins = v_cuts[v_argmins]
                                v_argmax = argmax(v_mins)
                                v_max = v_mins[v_argmax]
                            else
                                v_max = -1
                            end

                            # horizontal cut
                            if j > 1
                                h_cuts =  [[C[p,q,i,t], C[p,q+t,i,j-t]][k] for t=1:j-1, k=1:2]
                                for t=1:j-1
                                    if h_cuts[t,1]==inf || h_cuts[t,2]==inf
                                        h_cuts[t,:] = [-1, -1]
                                    end
                                end
                                h_argmins = argmin(h_cuts, dims = 2)
                                h_mins = h_cuts[h_argmins]
                                h_argmax = argmax(h_mins)
                                h_max = h_mins[h_argmax]
                            else
                                h_max = -1
                            end                    

                            mini = [D[p,q,i,j], v_max, h_max]

                            C[p,q,i,j] = maximum(mini)
                            argcuts[p,q,i,j] = argmax(mini)
                            if argcuts[p,q,i,j] == 2
                                idxcuts[p,q,i,j] = v_argmax[1]
                            elseif argcuts[p,q,i,j] == 3
                                idxcuts[p,q,i,j] = h_argmax[1]
                            end
                        end
                    end
                end
            end
        end
    end

    
    return C, P, D, midpoints_x, midpoints_y, argcuts, idxcuts
end