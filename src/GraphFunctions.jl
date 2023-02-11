function Base.IteratorSize(::Base.Iterators.Flatten{Vector{Combinatorics.MultiSetPermutations{Vector{Int64}}}})::Base.HasLength
    return Base.HasLength()
end
function Base.length(F::Base.Iterators.Flatten{Vector{Combinatorics.MultiSetPermutations{Vector{Int64}}}})::Int64
    return sum(length, F.it)
end

function get_weights(l::Int64, d::Int64)::Vector{Vector{Int64}}

    return collect(Iterators.Flatten([multiset_permutations(p, l) for p in partitions(d,l)]))
end
