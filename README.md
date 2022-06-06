# IntegratedOscillatorModel.jl
# Background
This module was made during my bachelors thesis. It implements a mathematical model developed 
the following publication [McKenna et. al 2016](https://doi.org/10.1016/j.bpj.2015.11.3526)
The specific model used here is the same as in [Marinelli et. al 2018](https://doi.org/10.1016/j.jtbi.2018.06.017)
The model describes the flux ot certain ions and molecules seen in the beta-cell in the human
pancreas. This is done by a non linear system of first order differential equations and it's
goal is to get insights about the insulin exostosis.

I will not explain the specific parameters of the model in this README.
Please refer to the above mentioned papers and the coments in the params.jl file.

If you want to gain a broader overview of beta-Cell modeling, I can recommend [this paper](https://doi.org/10.2337/dbi17-0004).

# Goals
The goal of this module is to enable easy manipulation of parameters within the model,
to obtain insides about it's viability and inner workings.
Fast calculation was also a concern so multithreading was implemented as well as
helper functions to deal with the arising problems.
A novel approach is, to include the simulation of active substances in the model.

# Usage
The main function to use, is the simulate function.
It accepts a System struct, that has the following structure:

## System
```julia
@with_kw struct System
    ode_func      :: Function = sys 
    y0            :: AbstractVecOrMat{Float64} = y0
    time          :: Integer = 60
    plot_args     :: Dict = Dict()
    params        :: Dict{String,Any} = deepcopy(params)
    plot_params   :: AbstractVector{Integer} = [1,2,3,4,5,6,7,8]
    change_params :: Dict{String,Number} = Dict()
    meds          :: Vector{Meds} = []
end
```
The above struct uses the [Parameters.jl](https://github.com/mauro3/Parameters.jl) package.

- ode_func\
    The function, that is used as the oder function. In most cases this parameter should stay the same.

- y0 \
    The starting values.
    Another set of values (called y0_stat) contains the values to achive a almost constatn solution
- time \
    The simulation time in minutes
- plot_args \
    Extra arguments, that are used in the automatic plot generation
- params \
    The parameterset of the simulation.
    This will be modified, so copy or deepcopy should be used to mittigate unwanted side effects.
- plot_params \
    The values of the ode system that shall be contained in the aut generated plot.
    1. V - Membrane potential
    2. N - Kalium
    3. Ca - Calcium
    4. Ca_er - Calcium (endoplasmatic reticulum)
    5. Ca_m - Calcium (mitochondria)
    6. ADP
    7. F6P
    8. FBP
- change_params \
    The key value pairs are used to modify some parameters given with the params field.
- meds \
    A list of Activa objects
## Activa
Activa represent active substances, that can be added to the system as further behavour modifiers.

```julia
@with_kw mutable struct Activa
    # starting time in minutes
    time :: Real
    # duration in inutes
    duration :: Number = Inf
    # dosage
    dose :: Float64
    name :: String = ""
    # will be used in en exponential function to simulate
    # exponential dosage decrease or decline
    fade :: Float64 = 1
    # (dose, current_ode_state) -> Float64
    # This function is used to calculate the effect of a active substance
    func :: Function = id
    # The name of the param, that is modified
    param :: String
end
```

## simulate
The interface to use, is the simulate function.
It accepts a Settings object and will run a simulation accordingly.
The return value is a tuple, containg an OdeSolution object, a matrix
containing the solution and a plot, of the solution.

```julia
using IntegratedOscillatorModel
ode_solution, solution_matrix, solution_plot = simulate(System())
```

It is also possible to simulate different modifications in one command.
In this scenario 3 dicts will be returned, containing the solutions, the matrices,
and the plots. In addition a plot, of all solutions will be returned as the fourth value
Because, this function uses the [@threads makro](https://docs.julialang.org/en/v1/manual/multi-threading/) to parralelize the operations, the SortedDict has a key for every value, given in the value array.

```julia
using IntegratedOscillatorModel
sol_dict, mat_dict, plot_dict, overview_plot = simulate(
    System(
        time = 120
        change_params = Dict(
            "vpdh" => 0.001
        )),
    "Jgk",
    [0.001, 0.05, 0.01]
    )
```

At last, two keys and two lists can be provided like in the following example:
```julia
using IntegratedOscillatorModel
sol_dict, mat_dict, plot_dict, overview_plot = simulate(
    System(
        time = 120
        change_params = Dict(
            "vpdh" => 0.001
        )),
    "Jgk",
    [0.001, 0.05, 0.01]
    "gca",
    [200, 400, 800]
    )
```

## simulate activa

Four types of activa are currently predefined:

| Function     | Active Substance | param |
| :---------- | :----------------- | :---- |
| Tolb         | Tolbuatmid | - gkatpbar |
| Dz           | Diazoxide | + gkatpbar |
| Glucose      | Glucose | + Jgk |
| KCl          | KalciumChloride | + vk |

#### Example

```julia
using IntegratedOscillatorModel
sol = simulate(System(
        time = 140,
        plot_params = [1,3,6],
        change_params = Dict(
            "vpdh" => 0.001,
            "Jgk" => 0.001,
            "gkca" => 800
        ),
        meds = [
            Tolb(40,  50,   duration = 20),
            Tolb(60,  100,  duration = 20),
            Tolb(80,  200,  duration = 20),
            Tolb(100, 500,  duration = 20),
            Tolb(120, 1000, duration = 20),
        ]
    ),
)
```
# Inner Workings and Hints

### Tolerance
Since the simulated ode system is quit fast oscillating, the sovler will assume an initial
tolerance of abstol = 1e-6 and reltol = 1e-6.
If the calculation failed because of a domain error (mostly a failed sqrt(-n)) the 
tolerances will be decreased by 6 orders of magniuted.
If the system can't by solved with a tolerance of 1e-18, the simulations failes and returns nothing.

### Medication implementation
A periodic callback is called every 100 steps in the simulation, and the medications are evaluated.
An evaluation at every step would be very slow and error prone.

### Reexports
The Plots and LaTeXStrings modules are reexported.

### Legacy
I originaly desinged the module to work with a named tuple of options.
This syntax is still accpeted, to achive compatibility with my bachelors thesis simulations.

