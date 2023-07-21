using Plots
using LazySets

function discrepancy(S::Matrix{<:Real},w::Vector{<:Integer},Π::Vector{<:Integer})
    k = maximum(Π)  
    SΠ_disc = [abs(sum((Π.==j) .* (w.==1)) - sum((Π.==j) .* (w.==-1))) for j=1:k]
    disc = minimum(SΠ_disc)
    return disc
end


function relabler!(Π::Array{<:Integer}, new_names::Vector{<:Integer})
    for i = 1:length(Π)
        Π[i] = new_names[Π[i]]
    end
    return Π
end


function partitions(n::Int, k::Int)
    if k == 1
        return [1 for i in 1:n]
    elseif n == k
        return [i for i in 1:n]
    else
        ## This is essentially the recurrence of Stirling numbers of the second kind. (https://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind)
        res = Array{Integer,2}(undef,n,0)
        for i in 1:k
            sub_par = relabler!(partitions(n-1, k), circshift([j for j in 1:k],i-1))
            par = vcat(ones(Int,1,size(sub_par,2)),sub_par)
            res = hcat(res,par)
        end
        sub_par = relabler!(partitions(n-1, k-1), [j for j in 2:k])
        par = vcat(ones(Int,1,size(sub_par,2)),sub_par)
        res = hcat(res,par)
        return res
    end
end


function coarseness_bf(S::Matrix{<:Real}, w::Vector{<:Integer})
    c = 1
    Π_opt = Array{Integer,2}(undef,n,0)

    for k = 1:n
        par = partitions(n,k)

        for p = 1:size(par,2)
            Π = par[:,p]

            hulls = Array{Array{Array{Float64, 1}, 1}, 1}(undef,k)
            for j = 1:k
                j_points = S[Π .== j,:]
                j_set = [j_points[i,:] for i in 1:size(j_points,1)]
                hulls[j] = convex_hull(j_set)
            end
            flag = true
            for j = 1:k-1
                for l = j+1:k
                    flag = flag * isempty( VPolygon(hulls[j]) ∩ VPolygon(hulls[l]) )
                end
            end

            if flag 
                d = discrepancy(S,w,Π)
                if d == c
                    Π_opt = hcat(Π_opt, Π)
                elseif d > c
                    Π_opt = Π
                    c = d
                end
            end
        end
    end
    
    return c, Π_opt
end