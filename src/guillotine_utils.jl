function sorted2midpoints(sorted)
    n = length(sorted)
    midpoints = zeros(n + 1)
    Δ = (sorted[1] + sorted[n])/n
    midpoints[1] = sorted[1] - Δ
    midpoints[n+1] = sorted[n] + Δ
    midpoints[2:n] = 0.5(sorted[2:n] + sorted[1:n-1])
    return midpoints
end