# default labels and starting values for simulation
const labels = [L"V" L"N" L"Ca" L"Ca_{er}" L"Ca_m" "ADP" "F6P" "FBP"];
const y0 = [-60; 0; 0.1; 185; 100; 780; 60; 40];
# y0_stat produces a stationary solution
const y0_stat = [ -60.486843436763024
                    0.0001367295835481462
                    0.06376844044216665
                  127.60067923018624
                   25.57114498828348
                  809.2663464858219
                    9.141692725854208
                    4.432631767186038e-6 ];


# defines a struct with default values
# A new system can be instancated with keyword arguments
#
# Example
# System(t = 120) -> Struct with default parameters but 120 min length
#
@with_kw struct System
    ode_func :: Function = sys  
    y0 :: AbstractVecOrMat{Float64} = y0
    time :: Integer = 60
    plot_args :: Dict = Dict()
    params :: Dict{String,Any} = deepcopy(params)
    plot_params :: AbstractVector{Integer} = [1,2,3,4,5,6,7,8]
    change_params :: Dict{String,Number} = Dict()
    meds :: Vector{Meds} = []
    logging :: Bool = true
    solver :: SolverType = Tsit5()
end


# This function attempts to solve the ode system.
# If an DomainError exception occures it retries with
# increased accuracy. If an unexpted exception is thrown
# it is printed to stdout
# Refere to the DifferentialEquations.jl docs
# https://diffeq.sciml.ai/stable/tutorials/ode_example/
#trysolve(system :: System, callback, iteration) = trysolve(to_named_tuple(trysolve), callback, iteration)
function trysolve(system, callback, iteration)
    tol = 1e-1^(iteration*6)
    backup = deepcopy(system) # Create copy for further iterations.
			      # This is necessary because the Med and ExpMed objects are destroyed 
    problem = ODEProblem(
        system.ode_func,
        system.y0,
        (0.0, (
            system.time < 1000)
                ? system.time * 6000
                : system.time),
        system.params)
    if tol <= 1e-19 # maximum tolerance
        return nothing
    end
    try
	    # the following callback is called every 100 timesteps and checks for medications to apply
        cb = (callback) ? PeriodicCallback(parse_medications, 100) : nothing
	    # https://diffeq.sciml.ai/stable/basics/common_solver_opts/
        return solve(problem,
                Tsit5(),
                saveat   = 1,
                maxiters = 1e6,
                reltol   = tol,
                abstol   = tol,
                callback = cb,
                progress = system.logging,
                progress_steps = 1,
        )
    catch err
        if isa(err, DomainError)
            return trysolve(backup, callback, iteration+1)
        else
            println(err)
        end
    end
end

function add_meds_to_plot(p, meds :: Vector{Union{Activa,ExpMed,Med}}, matrix)
    plot_max = matrix |> Iterators.flatten |> collect |> max
    p
    for (i, m) ∈ enumerate(filter(x -> x isa Activa, meds))
            plot!(
                [m.time, m.time + m.duration],
                [plot_max + 20 * i, plot_max + 20 * i],
                label = m.name,
                linestyle = :dash
                )
            annotate!(m.time + m.duration / 2, plot_max + 20 * (i + 0.4), string(m.dose))
            #println("med annotations cannot be added with legacy Med or ExpMed structs")
    end
    return p
end

function gen_solution_plot(solution :: ODESolution, system :: System)
    gr()
    # apply scaling factors to fit all parameters in one plot,
    # since they are numerically unreliable anyways
    matrix = scale_solution_columns(solution)
    p = plot(
        solution.t/6000, 
        matrix[system.plot_params];
        label = hcat(labels[system.plot_params]...),
        format=:png,
        fontfamily = "Computer Modern",
        titlefontsize  = 14,
        guidefontsize  = 14,
        tickfontsize   = 14,
        legendfontsize = 13,
        framestyle = :default,
        width=1.2,
        size = (900,600),
        palette = :tab10,
        grid = false,
        legend = (1.07,0.8),
        rightmargin = 42mm,
        system.plot_args...
        )
        
    print(p)
    p = add_meds_to_plot(p, system.meds, matrix)
    return p
end

# Takes a settings object (NamedTuple) and changes to parameters for the simulation
# accordingly. It returns DiffEq Solution object, the time series for the 8 values
# and a plot fo the solution.

function simulate(system :: System; iteration = 1)
    for (key,val) ∈ system.change_params
        system.params[key] = val;
    end
    system.params["meds"] = system.meds
    
    @time solution = trysolve(system, system.meds |> length > 0, iteration)
    
    (solution === nothing) ? nothing : (
        solution,
        rowtocolumn(solution.u),
        gen_solution_plot(solution,system)
    )
end

# maintina legacy compatablity with old named tuple args
simulate(system) = simulate(System(system...))

# can be called to put multiple vales for one parameter
# and produce one solution for each
# If multiple threads are available the simulations
# will be run in parallel

function loopvals(key::String, vals::Array, system)
    sols = Dict()
    sols_mat = Dict()
    plots = Dict()
    
    @threads for val ∈ vals
        local s = deepcopy(system)
        s.params[key] = val
        s.plot_args[:title] = key * " = " * string(val)
        try
            sol, mat, p = simulate(s)
            sols[val] = sol
            sols_mat[val] = mat
            plots[val] = p
        catch
            println("simulate " * key * " = " * string(val) * " failed")
        end
    end
    acc_plot = plot(
        sortedvalues(plots)...,
        dpi=300,
        label = hcat(labels[system.plot_params]...),
    )
    return (sols, sols_mat, plots, acc_plot)
end

# same as above but for two chaning variables
function loopvals(key1::String, vals1::Vector, key2::String, vals2::Vector, system)
    local sols::Dict{Any, Dict}     = Dict()
    local sols_mat::Dict{Any, Dict} = Dict()
    local plots::Dict{Any, Dict}    = Dict()
    
    @threads for val1 ∈ vals1
    #for val1 ∈ vals1
            sols[val1] = Dict()
            sols_mat[val1] = Dict()
            plots[val1] = Dict()
        @threads for val2 ∈ vals2
        #for val2 ∈ vals2
            local s = deepcopy(system)
            s.params[key1] = val1
            s.params[key2] = val2
            s.plot_args[:title] = key1 * " = " * string(val1) * ", " * key2 * " = " * string(val2)
            try
                sol, mat, p = simulate(s)
                sols[val1][val2] = sol
                sols_mat[val1][val2] = mat
                plots[val1][val2] = p
            catch
                println("simulate " * key1 * " = " * string(val1) * " " * key2 * " = " * string(val2) *" failed")
            end
        end
    end
    acc_plot = plot(
        sortedvalues(plots)...,
        dpi=300,
        label = hcat(labels[system.plot_params]...),
    )
    
    return (sols, sols_mat, plots, acc_plot)
end

# use multiple dispatch to avoid calls to loopvals
simulate(system :: System, key1::String, vals1::Vector, key2::String, vals2::Vector) = loopvals(key1, vals1, key2, vals2, system)
simulate(system :: System, key1::String, vals1::Vector) = loopvals(key1, vals1, system)
