# import functions explicitly to override them
import Base.sort
import Base.max
import Base.min

using DataStructures

# returns a list of n times the vec concatenated
function cycle(vec, n)
    return collect(Base.Iterators.flatten((Base.Iterators.repeated(vec, n))))
end

# calculate the nernst potential
function vkf(ke,ki)
    8.314 * 293.15/96485.332 * log(ke/ki) * 1000
end

function vkf(ki)
    8.314 * 293.15/96485.332 * log(15/ki+1) * 1000
end

# is used to get the get the intermediate stepts for variables used
# to calculate the changes in the main variables
function setupstats()
    return Dict([
        ("t", []),
        ("vk", []),
        ("gkatpbar", []),
        ("gca", []),
        ("gkca", []),
    ])
end

function savestats(i)
    global stats
    push!(stats["t"]    , i.t)
    push!(stats["vk"]    , i.p["vk"])
    push!(stats["gkatpbar"]    , i.p["gkatpbar"])
end

# takes a dictionary and sorts it alphanumericaly by keys (returns a SortedDict)
function dsort(d::Dict)
    sorted_keys = map(x -> x[1], sort(pairs(d) |> collect, by=x->x[1]))
    sd = SortedDict()
    map(x -> sd[x] = d[x], sorted_keys)
    return sd
end

function getkeys(v)
    try
        return keys(v)
    catch
        return fieldnames(typeof(v))
    end
end

# convert struct to named tuple
# https://discourse.julialang.org/t/get-fieldnames-and-values-of-struct-as-namedtuple/8991
gentup(struct_T) = NamedTuple{( fieldnames(struct_T)...,), Tuple{(fieldtype(struct_T,i) for i=1:fieldcount(struct_T))...}}
@generated function to_named_tuple(x)
    nt = Expr(:quote, gentup(x))
    tup = Expr(:tuple)
    for i=1:fieldcount(x)
        push!(tup.args, :(getfield(x, $i)) )
    end
    return :($nt($tup))
end

# calculates the difference between each element of a vector
function vecdiff(vec::Vector)
    map(v -> v[1] - v[2], zip(vec, [vec[2:end]; vec[end]]))
end

# enter the total time interval to get dy/dt value for the vector
function vecdiff(vec::Vector, interval::Number)
    vecdiff(vec) ./ (interval ./ length(vec))
end

# returns max/min of vector (there is also a builtin findmax/findmin function)
function max(vec::Vector)
    reduce(vec) do x,y
       x>y ? x : y end
end

function min(vec::Vector)
    reduce(vec) do x,y
       y>x ? x : y end
end

# returns a list of the values stored in a dict
function values(d::Dict; sorted=false)
    if sorted
        return sorted_values(d)
    end
    map(x -> x[2], collect(pairs(d)))
end
    
# returns a list of values in a dict sorted by key value
function sortedvalues(d::Dict)
    map(x -> x[2], sort(pairs(d) |> collect, by=x->x[1]))
end

function sortedvalues(d::Dict{Any, Dict})
    reduce(vcat, map(x -> sortedvalues(x), map(x -> x[2], sort(pairs(d) |> collect, by=x->x[1]))))
end

# By default the numerical integration returns an array for each timestep
# that contains the values for all variables at that time.
# For plotting and further calculations its easier to have a seperate
# time series for each variable.
function rowtocolumn(vec::AbstractVector{T}) where T <: AbstractVector
    dim1 = length(vec)
    dim2 = length(vec[1])
    vects = Vector{Vector}(undef, dim2)
    map(x -> vects[x] = Array{eltype(vec[1])}(undef, dim1), 1:dim2)
    @inbounds @fastmath for i in 1:dim1, j in 1:dim2
        vects[j][i] = vec[i][j]
    end
    return vects
end

# obsolete, but could be usefull for some calculations
# this converts a vector of vectors to a matrix.
function vector_to_matrix(vecvec::AbstractVector{T}) where T <: AbstractVector
    dim1 = length(vecvec)
    dim2 = length(vecvec[1])
    my_array = Array{eltype(vecvec[1]), 2}(undef, dim1, dim2)
    @inbounds @fastmath for i in 1:dim1, j in 1:dim2
        my_array[i,j] = vecvec[i][j]
    end
    return my_array
end

function vec_to_matrix(vecvec::AbstractVector{T}) where T <: AbstractVector
    return vector_to_matrix(vecvec)
end

# scales the solutions to all be visible in the same plot
function scale_solution_columns(sol)
    #mat = vec_to_matrix(sol.u)
    mat = rowtocolumn(sol.u)
    mat[2] *= 100
    mat[3] *= 750
    mat[6] /= 7
    #if size(mat)[1] == 9
        #mat[9] *= 10
    #end
    return mat
end
