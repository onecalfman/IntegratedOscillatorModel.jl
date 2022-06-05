__precompile__()

module IntegratedOscillatorModel

print("1.0.3")

# import functions explicitly to override them
import Base.sort
import Base.max
import Base.min

# library imports
using Base.Threads
using DataStructures
using Reexport
using DifferentialEquations
using Parameters
@reexport using LaTeXStrings
@reexport using Plots; gr();
@reexport using Plots.Measures

include("plot_settings.jl")
include("helper.jl")
include("meds.jl")
include("params.jl")
include("ode_systems.jl")
include("solver_utils.jl")

# exported functions and variables that are available after using IntegratedOscillatorModel
export simulate
export loopvals
export labels
export y0
export y0_stat
export sys
#export pfk_activity
export params
export sortedvalues
export rowtocolumn
export max
export min
#export vecdiff
export vkf
#export cycle
export Med
export ExpMed
export Meds
export trysolve
export sys_const_pfk
export System

end # module
