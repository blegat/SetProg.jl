function SumOfSquares.moment_matrix(cref::SetConstraintRef)
    return SumOfSquares.moment_matrix(transformed_constraint(cref))
end
