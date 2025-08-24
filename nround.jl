"""
  nround(x)

This function returns the correctly-rounded integer value (of `Int` type)
for a given value of x, which should be of `Real` type.
"""
function nround(x)
# This function returns the correct rounded integer value for a given value of x.
# Sanity check first.
if typeof(x) <: Real
  return Float64(x) - floor(Float64(x)) < 0.5 ? Int(floor(Float64(x))) : Int(ceil(Float64(x)))
else
  error("Cannot use this function because the input argument is of type ", typeof(x))
end
end
