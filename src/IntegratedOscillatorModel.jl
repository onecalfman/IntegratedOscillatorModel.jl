__precompile__()

module IntegratedOscillatorModel
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
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

@reexport using LaTeXStrings
@reexport using Plots
@reexport using Plots.Measures

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
export y0_analyzer_
export y0_stat
export sys
export sys_analyzer
#export pfk_activity
export df_from_sol
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

export Activa
export Dz
export Tolb
export Tg
export KCl
export Glucose

export Meds
export trysolve
export sys_const_pfk
export System

end # module
