using Plots
using LazySets

function plot_rb_points(
        S::Matrix{<:Real}, 
        w::Vector{<:Integer}, 
        Π::Vector{<:Integer}=[0]
    )
    
    S = Float64.(S)
    red_points = S[w.==1, :]
    blue_points = S[w.==-1, :]

    max_x, max_y = maximum(S, dims = 1)
    min_x, min_y = minimum(S, dims = 1)
    

    
    max_x = maximum(S[:,1])
    plt = scatter(size=(500, 500), legend=:none, 
                  xlims=(min_x-1, max_x+1),
                  ylims=(min_y-1, max_y+1), 
                  axis=0, framestyle=:box,
                  grid=0, color=:white)
    
    scatter!(plt, red_points[:, 1], red_points[:, 2], color=:red)
    scatter!(plt, blue_points[:, 1], blue_points[:, 2], color=:blue)
    
    
    
    if Π != [0]
        k = maximum(Π)        
        for j = 1:k
            j_points = S[Π .== j,:]
            j_set = [j_points[i,:] for i in 1:size(j_points,1)]
            j_hull = convex_hull(j_set)
            plot!(VPolygon(j_hull), color=:gray, alpha=0.2)
        end
        SΠ_disc = [abs(sum((Π.==j) .* (w.==1)) - sum((Π.==j) .* (w.==-1))) for j=1:k]
        
        disc = minimum(SΠ_disc)
        annotate!((.25, .9), text("Discrepancy: min $SΠ_disc = $disc", 8))
    end
    
    return plt
end
