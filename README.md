# IntegratedOscillatorModel.jl
# Background
This module was made during my bachelors thesis. It implements a mathematical model developed 
the following publication [McKenna et. al 2016](https://doi.org/10.1016/j.bpj.2015.11.3526)
The specific model used here is the same as in [Marinelli et. al 2018](https://doi.org/10.1016/j.jtbi.2018.06.017)
The model describes the flux ot certain ions and molecules seen in the beta-cell in the human
pancreas. This is done by a non linear system of first order differential equations and it's
goal is to get insights about the insulin exostosis.

# Goals
The goal of this module is to enable easy manipulation of parameters within the model,
to obtain insides about it's viability and inner workings.
Fast calculation was also a concern so multithreading was implemented as well as
helper functions to teal with the arising problems.

# How to use it
## Installation
This module can be installed from github.
Either
1. Enter the package mode by pressing ] in the julia repl
2. type:

        add https://github.com/onecalfman/IntegratedOscillatorModel.jl
or

    using Pkg
    Pkg.add(url="https://github.com/onecalfman/IntegratedOscillatorModel.jl")
## Usage

To use the module you have to create a settings object like in this example.

    settings = (
        ode_func      = sys,
        y0            = y0,
        time          = 140 * 6000,
        plot_args     = Dict(),
        params        = deepcopy(params),
        plot_params   = [1,3,6],
        change_params = Dict(
               "vpdh" => 0.001,
               "Jgk" => 0.001,
               "gkca" => 800,
            ),
        meds = [
            Med(20 + 20, "tolbutamid",  50),
            Med(20 + 40, "tolbutamid",  100),
            Med(20 + 60, "tolbutamid",  200),
            Med(20 + 80, "tolbutamid",  500),
            Med(20 + 100, "tolbutamid",  1000),
            ]
    );
    
    s, m, p = simulate(settings)

- ode_func
    This parameter is the function to numerical integrate.
    One could define variations of the here implemented sys function to create new behavior.
- y0
    starting conditions for the ode system.
    y0 is also a variable with the default parameters.
    In addition a y0_stat vector exists which has the initial settings to stay in a steady state
- time
    time in 100ms intervals
- plot_args
    currently defunct (should in theory allow to pass arguments for plotting)
- params
    the params to use for the calculation.
    Note: it is important to deepcopy the default params since julia passes objects by reference
    and the can be changed or deleted if passed directly
- plot_params
        a list of number representing the quantities that shall be plotted.
        the numbers represent (in this order)
        1. membrane potential
        2. kalium
        3. intracellular calcium
        4. calcium (ectoplasmatic reticulum)
        5. calcium (mitochondria)
        6. adp
        7. f6p
        8. fbp
- change_params
    a dict which overwrites the given parameters temporarily
- meds
    accepts list of Med or ExpMed structs.
    
## Meds
To simulate different drugs and their effect on the beta cell, Med or ExpMed objects can be passed to the system.
The consist of 

    struct Med            
        time::Real        # time in minutes at which the medication is added
        type::String      # type of medication (if type is equal to a key in the params dict that value will be changed to the value of dose)
        dose::Float64     # dose of medication (varies with type -> refer to handle_medication function) 
        fade::Float64     # time * fade^(current_time - time) -> sets decline (fade only exists on ExpMed)
    end

The available stimulants are
- tolbutamid
- gkatpbar (maximal conductance of atp dependent k channels)
- Dz (Diazoxide)
- glucose
- kserca (serca pump rate ca2+ er)
- KCl
- nifedipin (needs better fine tuning)
- any key contained in the params dict

Additional types of stimulants can be added in the "handle_medication" function.

## Simulate
The simulate function takes a settings dict or named tuple.
It returns a tuple of the form

    (ODESolution, matrix of solutions, plot of solution)
    
The simulate function tries to solve first with lower accuracy to increase speed.
If it fails it will try to solve with higher accuracy.
If the ode system can't be solved with an accuracy value of 1e-19
the solver fails and return nothing.

## Loopvals
To make variation of parameters easier the loopvals function takes in a string with the parameter
name to change and a vector of values that should be plugged in.
It can also accept two strings and two vectors. It then tries all possible combinations of these parameters.
The fourth object returned here is a basic plot containing all subplots.

    Jgk = collect(0.0003:0.0005:0.005)
    sols, mats, plots, _ = loopvals("Jgk", Jgk, settings)
    
    # or

    gkca = [50, 800]
    vpdh = [0.002, 0.009]
    sols, mats, plots, p = loopvals("gkca", gkca, "vpdh", vpdh, settings)

