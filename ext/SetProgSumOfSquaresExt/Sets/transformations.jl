
function zero_eliminate(set::Householder{T}, I) where {T}
    J = setdiff(1:dimension(set), I)
    p = MP.subs(set.p, set.x[I] => zeros(T, length(I)))
    return Householder(UnknownSet{T}(), p, set.h[J], set.z, set.x[J])
end
