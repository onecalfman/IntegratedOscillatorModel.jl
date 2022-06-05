# import functions explicitly to override them
import Base.sort
import Base.max
import Base.min

# library imports
using Base.Threads
using LaTeXStrings
using DifferentialEquations
using Plots; gr();
using Plots.Measures
using DataStructures
