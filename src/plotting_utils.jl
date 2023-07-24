using Plots
using LazySets

function plot_rb_points(
        S::Matrix{<:Real}, 
        w::Vector{<:Integer}, 
        Π::Vector{<:Integer}=[0];
        ann_disc = false,
        fig_size = (500, 500)
    )
    
    S = Float64.(S)
    red_points = S[w.==1, :]
    blue_points = S[w.==-1, :]

    max_x, max_y = maximum(S, dims = 1)
    min_x, min_y = minimum(S, dims = 1)
    
    range_x = max_x - min_x
    range_y = max_y - min_y
    
    plt = scatter(size=fig_size, legend=:none, 
                  xlims=(min_x-.2range_x, max_x+.2range_x),
                  ylims=(min_y-.2range_y, max_y+.2range_y), 
                  axis=0, 
                  framestyle=:none,
                  grid=0,
                  color=:white)
    
    MS = 6
    scatter!(plt, red_points[:, 1], red_points[:, 2], color=:red, markersize = MS)
    scatter!(plt, blue_points[:, 1], blue_points[:, 2], color=:blue, markersize = MS)
    
    
    
    if Π != [0]
        k = maximum(Π)        
        for j = 1:k
            j_points = S[Π .== j,:]
            j_set = [j_points[i,:] for i in 1:size(j_points,1)]
            j_hull = convex_hull(j_set)
            plot!(VPolygon(j_hull), color=:gray, alpha=0.2)
        end
        if ann_disc
            SΠ_disc = [abs(sum((Π.==j) .* (w.==1)) - sum((Π.==j) .* (w.==-1))) for j=1:k]

            disc = minimum(SΠ_disc)
            annotate!((.25, .9), text("Discrepancy: min $SΠ_disc = $disc", 8))
        end
    end
    
    return plt
end

################################################################################################

function plot_add_lines!(plt, S)
    function sorted2midpoints(sorted)
        n = length(sorted)
        midpoints = zeros(n + 1)
        Δ = (sorted[1] + sorted[n])/(n+1)
        midpoints[1] = sorted[1] - Δ
        midpoints[n+1] = sorted[n] + Δ
        midpoints[2:n] = 0.5(sorted[2:n] + sorted[1:n-1])
        return midpoints
    end
    
    sortperm_x = sortperm(S[:, 1])
    sorted_x = S[sortperm_x, 1]
    midpoints_x = sorted2midpoints(unique(sorted_x))

    sortperm_y = sortperm(S[:, 2])
    sorted_y = S[sortperm_y, 2]
    midpoints_y = sorted2midpoints(unique(sorted_y))
    
    vline!(midpoints_x, color=:gray, alpha = .3)
    hline!(midpoints_y, color=:gray, alpha = .3)
    
    return plt
end

################################################################################################

function plot_guillotine_line!(plt,p,q,i,j,k=0)
    k_max = 12
    try 
    rec_k[1] = max(rec_k[1],k)
    catch
    end
    k = k + 1
    if k < k_max
        if argcuts[p,q,i,j] == 2
            s = idxcuts[p,q,i,j]
            plot!(
                [midpoints_x[p+s], midpoints_x[p+s]],
                [midpoints_y[q],   midpoints_y[q+j]], 
                color = :black
            	)
            scatter!(
                [midpoints_x[p+s]],
                [.5midpoints_y[q]+.5midpoints_y[q+j]],
                markersize = 6,
                markerstrokewidth = 0,
                color = :white,
                alpha = .6
            )
            annotate!(
                midpoints_x[p+s],
                .5midpoints_y[q]+.5midpoints_y[q+j],
                text("$k",6,:black)
            )
            plot_guillotine_line!(plt,p,q,s,j,k)
            plot_guillotine_line!(plt,p+s,q,i-s,j,k)
        elseif argcuts[p,q,i,j] == 3
            t = idxcuts[p,q,i,j]
            plot!(
                [midpoints_x[p],   midpoints_x[p+i]],
                [midpoints_y[q+t], midpoints_y[q+t]], 
                color = :black
            )
            scatter!(
                [.5midpoints_x[p]+.5midpoints_x[p+i]],
                [midpoints_y[q+t]],
                markersize = 6,
                markerstrokewidth = 0,
                color = :white,
                alpha = .6
            )
            annotate!(
                .5midpoints_x[p]+.5midpoints_x[p+i],
                midpoints_y[q+t],
                text("$k",6,:black)
            )
            plot_guillotine_line!(plt,p,q,i,t,k)
            plot_guillotine_line!(plt,p,q+t,i,j-t,k)
        end
    end
end

################################################################################################
# utils for debugging
################################################################################################

function print_C_pyramid(C,label="",io="")
    n = size(C)[1]
    pad = n>9 ? 5 : 4
    println(io,"\n"*label)
    for i = 1:n
        println(io,lpad("i=$i ",pad+1," "),"-"^((n+2)*pad))
        for j = 1:n
            println(io,lpad("j=$j ",pad+1," ")," "^((n+1)*pad),"□$i×$j")
            for q = n-j+1:-1:1
                print(io,lpad("q=$q ",2*pad," "))
                for p = 1:n-i+1
                    print(io,lpad(C[p,q,i,j], pad, " "))
                end
                println(io,"")
            end
            println(io," "^(2*pad),[lpad("p=$p",pad," ") for p = 1:n-i+1]...)
        end
        println(io,"")
    end
    println(io,lpad("end ",pad+1," "),"-"^((n+2)*pad))
end

################################################################################################

function print_S_info(S,p,q,i,j)
    sortperm_x = sortperm(S[:, 1])
    sorted_x = S[sortperm_x, 1]
    midpoints_x = sorted2midpoints(sorted_x)

    sortperm_y = sortperm(S[:, 2])
    sorted_y = S[sortperm_y, 2]
    midpoints_y = sorted2midpoints(sorted_y)
    
    M = Matrix{Any}(undef,n,4); M[:,1] = 1:n; M[:,2:3]= round.(S,digits=2); foo(x) = x==1 ? 'r' : 'b'; M[:,4] = foo.(w) 
    display(M)

    println("\np = $p, q = $q, i = $i, j = $j\n--------------------------")
    print("                      "); for i=1:n+1 print(i,"  ") end; print("  midpoint index\n")
    println("sortperm_x          = $sortperm_x")
    println("sortperm_x[p:p+i-1] = $(sortperm_x[p:p+i-1])\n")
    print("                      "); for i=1:n+1 print(i,"  ") end; print("  midpoint index\n")
    println("sortperm_y          = $sortperm_y")
    println("sortperm_y[q:q+j-1] = $(sortperm_y[q:q+j-1])\n")


    int = intersect(sortperm_x[p:p+i-1],sortperm_y[q:q+j-1])
    println("intersection = $int\n")

    print("     "); for i=1:n print((w[i]==-1 ? " " : ""),i,"  ") end; print("  point index\n")
    println("w = $w\n" )
    println("Disc_pqij  = abs(sum(w[intersection])) = $(abs(sum(w[int])))")
    println("Count_pqij = length(intersection)      = $(length(int))")

    plt = plot_rb_points(S, w)
    plot_add_lines!(plt, S)
    xx = [midpoints_x[p], midpoints_x[p], midpoints_x[p+i], midpoints_x[p+i]]
    yy = [midpoints_y[q], midpoints_y[q+j], midpoints_y[q+j], midpoints_y[q]]
    plot!(Shape(xx,yy), color=:gray, alpha = .2)
    return plt
end