__precompile__()

module IntegratedOscillatorModel

# exported functions and variables that are available after using IntegratedOscillatorModel
export simulate
export loopvals
export labels
export y0
export y0_stat
export sys
export pfk_activity
export params
export sortedvalues
export rowtocolumn
export max
export min
export vecdiff
export vkf
export cycle
export Med
export ExpMed
export trysolve
export sys_const_pfk
export Plots
export plot
export mm


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

# settings default plot parameters
# this doesn't change anything and ist just kept for completnes
Plots.theme(:default)
default(dpi = 300,
 format=:png,
 fontfamily = "Computer Modern",
 titlefontsize  = 14,
 guidefontsize  = 14,
 tickfontsize   = 14,
 legendfontsize = 13,
 framestyle = :default,
 width=1.2,
 palette = :tab10,
 grid = false,
 legend=false
)

const PLOTS_DEFAULTS = Dict(
    :theme=>:default,
    :fontfamily => "Computer Modern")

# The Med and Expmed are used to specify medication usage in the simulation
struct Med            
    time::Real        # time in minutes at wich the medication is added
    type::String      # type of medication (if type is equal to a key in the params dict that value will be changed to the value of dose)
    dose::Float64     # dose of medication (varies with type -> refere to handle_medication function) 
end
    
mutable struct ExpMed # as above but with exponential decline
    time::Real
    type::String
    dose::Float64
    fade::Float64     # time * fade^(current_time - time) -> sets decline
end

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

#default params for the simulations
const params = Dict(
    "gca"      => 1000,     # maximal conductance ca membran channel
    "gkca"     => 50,       # maximal conducatnce ca dependet k channels
    "gk"       => 2700,     # maxiaml conductance of k ? 
    "gkatpbar" => 25000,    # maximal conductance of atp dependent k channels
        
    "vca"      => 25,       # nernst potential ca2+
    "vm"       => -20,      # shapre parameter ca2+ activation
    "sm"       => 12,       # shape parameter ca2+ activation
        
    "vk"       => -75,      # nernst potential k
    "k_ext"    => 30,       # extracellular k
    "k_i_ext"  => 30,       # initial extracellular k
    "k_int"    => 280,      # intracellular k
        
    "sn"       => 5,        # shape parameter k
    "taun"     => 20,       # time constant rectifying current
    "nin"      => -16,      # shape parameter rectifying current
    "kd"       => 0.5,      # MgADP^2 factor KATP activation function 
    "kdd"      => 17,       # MgADP factor KATP activation function
    "ktd"      => 26,       # ADP3- factor KATP activation function
    "ktt"      => 1,        # ATP4- factor KATP activation function
    "atot"     => 3000,     # apt + adp
    "alpha"    => 5.18e-18, # current to ion flux conversion factor
    "vcyt"     => 1.15e-12, # volume cytosol
    "kpmca"    => 0.2,      # pmca pump rate ca2+ cell membrane
    "kserca"   => 0.4,      # serca pump rate ca2+ er
    "pleak"    => 2.0e-4,   # leak across er membrane
    "kuni"     => 0.4,      # rate of mitochonria uniporters
    "knaca"    => 0.001,    # activity of mitochndrial ca na exchanger NCLX
    "Cm"       => 5300,     # membrane capacitance
    "fca"      => 0.01,     # fraction of ca2+ not bound to buffers
    "sigmam"   => 100,      # ratio of mitochonidria to cytosol volume
    "sigmaer"  => 31,       # ratio of ER to cytosol volume
    "Jgk"      => 0.001,    # glucokinase flux
        
    "k1"       => 30,       # params for pfk activity
    "k2"       => 1,        # params for pfk activity
    "k3"       => 50000,    # params for pfk activity 
    "k4"       => 1000,     # params for pfk activity
    "famp"     => 0.02,     # params for pfk activity
    "fatp"     => 20,       # params for pfk activity
    "ffbp"     => 0.2,      # params for pfk activity 
    "fbt"      => 20,       # params for pfk activity
    "fmt"      => 20,       # params for pfk activity
    "kpfk"     => 0.06,     # params for pfk activity
        
    "vpfk"     => 0.01,     # flux rate of pfk
    "taua"     => 300000,   # d ADP/dt time constant
    "kCaPDH"   => 200,      # factor for pdh activation function
    "vpdh"     => 0.009,    # flux rate of pdh
    
    "dummy"     => []       # set type to Dict{String, Any} (in julia 1.6.3 setting the global variable type
                            # is not supported. Therefore i use this hack, so i can add a Vector as a dict
                            # entry later on.
    );

# calculates the effects of medications
# if medication behavour shall be different 
# this is the place to change it/ add a new medication
function handle_medication(med, i::Any) 
    if med.type ≡ "tolbutamid"
        if med.dose > 0
                i.p["gkatpbar"] = params["gkatpbar"] - 0.02 * 1250 * med.dose
            else
                i.p["gkatpbar"] = params["gkatpbar"]
            end
    elseif med.type ≡ "Dz"
        if med.dose > 0
            i.p["gkatpbar"] = params["gkatpbar"] + 0.02 *  1250 * med.dose
        else
            i.p["gkatpbar"] = params["gkatpbar"]
        end
    elseif med.type ≡ "Tg"
        if med.dose > 0
                i.p["kserca"] = params["kserca"] - 0.04 * med.dose
            else
                i.p["kserca"] = params["kserca"]
            end
    elseif med.type ≡ "KCl"
            i.p["vk"] = vkf((med.dose > 4.8) ? med.dose : 4.8, 130)
    elseif med.type ≡ "nifedipin"
        Bool(med.dose) ? i.p["gca"] = 500 : i.p["gca"] = 1000
    elseif med.type ≡ "glucose"
        i.p["Jgk"] = med.dose / 1000
    elseif med.type ⊆ keys(params)
        i.p[med.type] = med.dose
    end
end

# wrapper functions for handle_medication wich determain
function calc_dosage(m::Med, i)
    handle_medication(m, i)
    return true
end

# this version determains if the expoentialy falling
# dosage is still over the minimal threshold
function calc_dosage(m::ExpMed, i)
    m.dose *= m.fade^((i.t - m.time*6000)/6000)
    #println(m.dose)
    handle_medication(m, i)
    if m.dose ≤ 0.01
        m.dose = 0
        return true
    elseif m.type == "KCl" && m.dose ≤ 4.8
        return true
    end
    return false
end

function parse_medications(i)
    try
        savestats(i)
    catch
        print("")
    end
    for m ∈ i.p["meds"]
        if 6000 * m.time ≤ i.t
            if calc_dosage(m, i);
                filter!(n -> n ≠ m, i.p["meds"])
            end
        end
    end
end

# pfk activity function 
function pfk_activity(atp,adp,f6p,fbp,amp, params)
	# (alpha,beta,gamma,delta);
	# (0,0,0,0);
    
	weight1=1;
	topa1=0;
	bottom1=1;
	
	# (0,0,0,1);
	weight2=atp^2/params["k4"];
	topa2=topa1;
	bottom2=bottom1+weight2;
	
	# (0,0,1,0);
	weight3=f6p^2/params["k3"];
	topa3=topa2+weight3;
	bottom3=bottom2+weight3;
	
	# (0,0,1,1);
	weight4=(f6p*atp)^2/(params["fatp"]*params["k3"]*params["k4"]);
	topa4=topa3+weight4;
	bottom4=bottom3+weight4;
	
	# (0,1,0,0);
	weight5=fbp/params["k2"];
	topa5=topa4;
	bottom5=bottom4+weight5;
	
	# (0,1,0,1)
	weight6=(fbp*atp^2)/(params["k2"]*params["k4"]*params["fbt"]);
	topa6=topa5;
	bottom6=bottom5+weight6;
	
	# (0,1,1,0)
	weight7=(fbp*f6p^2)/(params["k2"]*params["k3"]*params["ffbp"]);
	topa7=topa6+weight7;
	bottom7=bottom6+weight7;
	
	# (0,1,1,1)
	weight8=(fbp*f6p^2*atp^2)/(params["k2"]*params["k3"]*params["k4"]*params["ffbp"]*params["fbt"]*params["fatp"]);
	topa8=topa7+weight8;
	bottom8=bottom7+weight8;
	
	# (1,0,0,0);
	weight9=amp/params["k1"];
	topa9=topa8;
	bottom9=bottom8+weight9;
	
	# (1,0,0,1);
	weight10=(amp*atp^2)/(params["k1"]*params["k4"]*params["fmt"]);
	topa10=topa9;
	bottom10=bottom9+weight10;
	
	# (1,0,1,0);
	weight11=(amp*f6p^2)/(params["k1"]*params["k3"]*params["famp"]);
	topa11=topa10+weight11;
	bottom11=bottom10+weight11;
	
	# (1,0,1,1);
	weight12=(amp*f6p^2*atp^2)/(params["k1"]*params["k3"]*params["k4"]*params["famp"]*params["fmt"]*params["fatp"]);
	topa12=topa11+weight12;
	bottom12=bottom11+weight12;
	
	# (1,1,0,0)
	weight13=(amp*fbp)/(params["k1"]*params["k2"]);
	topa13=topa12;
	bottom13=bottom12+weight13;
	
	# (1,1,0,1);
	weight14=(amp*fbp*atp^2)/(params["k1"]*params["k2"]*params["k4"]*params["fbt"]*params["fmt"]);
	topa14=topa13;
	bottom14=bottom13+weight14;
	
	# (1,1,1,0) -- the most active state of the enzyme;
	weight15=(amp*fbp*f6p^2)/(params["k1"]*params["k2"]*params["k3"]*params["ffbp"]*params["famp"]);
	topa15=topa14;
	topb=weight15;
	bottom15=bottom14+weight15;
	
	# (1,1,1,1);
	weight16=(amp*fbp*f6p^2*atp^2)/(params["k1"]*params["k2"]*params["k3"]*params["k4"]*params["ffbp"]*params["famp"]*params["fbt"]*params["fmt"]*params["fatp"]);
	topa16=topa15+weight16;
	bottom16=bottom15+weight16;
	
	Jpfk= params["vpfk"]*(topb + params["kpfk"]*topa16)/bottom16;
end

# the main ode function
function sys(dy, y, params, t)
    #v, n, c, cer, cam, adp, f6p, fbp = y

    # calc atp/adp ratios
    rad     = sqrt(-4*y[6]^2+(params["atot"]-y[6])^2);
    atp     = (params["atot"]+rad-y[6])/2;
    mgadp   = 0.165*y[6]; # magnesium adp complex
    adp3m   = 0.135*y[6]; # adp for transport into matrix  (löffler petrides s. 238)
    atp4m   = 0.05 *atp;  # atp for transport into matrix (löffler petrides s. 238)
    amp     = y[6]^2/atp; # adenosine monophosphate
    
    # flux activation functions
    topo    = 0.08+0.89*mgadp^2/params["kdd"]^2+0.16*mgadp/params["kdd"] ;
    bottomo = (1+mgadp/params["kdd"])^2*(1+atp4m/params["ktt"] + adp3m/params["ktd"]) ;
    katpo   = topo/bottomo;                                                 # I_K(ATP) activation function
    minf    = 1/(1+exp((params["vm"]-y[1])/params["sm"]));                  # ca pmca activation function
    ninf    = 1/(1+exp((params["nin"]-y[1])/params["sn"]));                 # rectifing current acti
    qinf    = y[3]^2/(params["kd"]^2+y[3]^2);                               # I_K(Ca) activation functin
    # fluxes                                                                   
    
    ik      = params["gk"]*y[2]*(y[1]-params["vk"]);                        # rectifying current
    ikca    = -params["gkca"]*qinf*(params["vk"]-y[1]);                     # ca dependent k current
    ikatp   = params["gkatpbar"]*katpo*(y[1]-params["vk"]);                 # k flux atp dependent
    ica     = params["gca"]*minf*(y[1]-params["vca"]);                      # ca flux for action potentials
    
    Jer     = params["kserca"]*y[3] - params["pleak"]*(y[4]-y[3]);          # ca flux density across er membrane
    Jm      = params["kuni"]*y[3] - params["knaca"]*(y[5]-y[3]);            # ca flux density across mitochonidium
    Jmem    = -(params["alpha"]/params["vcyt"]*ica + params["kpmca"]*y[3]); # ca flux density (cell membrane)
                                                                               
    # glycolitic oscillations                                                                                       
    Jpfk    = pfk_activity(atp, y[6], y[7], y[8], amp, params)              # pfk activity
    sinfty  = y[5]/(y[5]+params["kCaPDH"]);                                 # Michaelis-Menten function   
    Jpdh    = params["vpdh"]*sinfty*sqrt(y[8]);                             # pdh activity
    
    #save_stats(amp, atp, katpo, minf, ninf, qinf, ik, ikca, Jer, ikatp, Jm, ica, Jmem, Jpfk, sinfty, Jpdh);
    dy[1] = -(ica + ik + ikca + ikatp)/params["Cm"]
    dy[2] = -(y[2]-ninf)/params["taun"]
    dy[3] = params["fca"]*(Jmem - Jm - Jer)
    dy[4] = params["fca"]*params["sigmaer"]*Jer
    dy[5] = params["fca"]*params["sigmam"]*Jm
    dy[6] = (atp-exp((1+2.2 * Jpdh/(Jpdh+0.05)) * (1-y[3]/0.35))*y[6])/params["taua"]
    dy[7] = 0.3*(params["Jgk"]-Jpfk)
    dy[8] = Jpfk-Jpdh/2
end

# the main ode function with constant pfk
function sys_const_pfk(dy, y, params, t)
    #v, n, c, cer, cam, adp, f6p, fbp = y

    # calc atp/adp ratios
    rad     = sqrt(-4*y[6]^2+(params["atot"]-y[6])^2);
    atp     = (params["atot"]+rad-y[6])/2;
    mgadp   = 0.165*y[6]; # magnesium adp complex
    adp3m   = 0.135*y[6]; # adp for transport into matrix  (löffler petrides s. 238)
    atp4m   = 0.05 *atp;  # atp for transport into matrix (löffler petrides s. 238)
    amp     = y[6]^2/atp; # adenosine monophosphate
    
    # flux activation functions
    topo    = 0.08+0.89*mgadp^2/params["kdd"]^2+0.16*mgadp/params["kdd"] ;
    bottomo = (1+mgadp/params["kdd"])^2*(1+atp4m/params["ktt"] + adp3m/params["ktd"]) ;
    katpo   = topo/bottomo;                                                 # I_K(ATP) activation function
    minf    = 1/(1+exp((params["vm"]-y[1])/params["sm"]));                  # ca pmca activation function
    ninf    = 1/(1+exp((params["nin"]-y[1])/params["sn"]));                 # rectifing current acti
    qinf    = y[3]^2/(params["kd"]^2+y[3]^2);                               # I_K(Ca) activation functin
    # fluxes                                                                   
    
    ik      = params["gk"]*y[2]*(y[1]-params["vk"]);                        # rectifying current
    ikca    = -params["gkca"]*qinf*(params["vk"]-y[1]);                     # ca dependent k current
    ikatp   = params["gkatpbar"]*katpo*(y[1]-params["vk"]);                 # k flux atp dependent
    ica     = params["gca"]*minf*(y[1]-params["vca"]);                      # ca flux for action potentials
    
    Jer     = params["kserca"]*y[3] - params["pleak"]*(y[4]-y[3]);          # ca flux density across er membrane
    Jm      = params["kuni"]*y[3] - params["knaca"]*(y[5]-y[3]);            # ca flux density across mitochonidium
    Jmem    = -(params["alpha"]/params["vcyt"]*ica + params["kpmca"]*y[3]); # ca flux density (cell membrane)
                                                                               
    # glycolitic oscillations                                                                                       
    #Jpfk    = pfk_activity(atp, y[6], y[7], y[8], amp, params)              # pfk activity
    Jpfk    = 0.003
    sinfty  = y[5]/(y[5]+params["kCaPDH"]);                                 # Michaelis-Menten function   
    Jpdh    = params["vpdh"]*sinfty*sqrt(y[8]);                             # pdh activity
    
    #save_stats(amp, atp, katpo, minf, ninf, qinf, ik, ikca, Jer, ikatp, Jm, ica, Jmem, Jpfk, sinfty, Jpdh);
    dy[1] = -(ica + ik + ikca + ikatp)/params["Cm"]
    dy[2] = -(y[2]-ninf)/params["taun"]
    dy[3] = params["fca"]*(Jmem - Jm - Jer)
    dy[4] = params["fca"]*params["sigmaer"]*Jer
    dy[5] = params["fca"]*params["sigmam"]*Jm
    dy[6] = (atp-exp((1+2.2 * Jpdh/(Jpdh+0.05)) * (1-y[3]/0.35))*y[6])/params["taua"]
    dy[7] = 0.3*(params["Jgk"]-Jpfk)
    dy[8] = Jpfk-Jpdh/2
end

# the main ode function with fixed value for katpo variable
function sys_katpo(dy, y, params, t)
    #v, n, c, cer, cam, adp, f6p, fbp = y

    # calc atp/adp ratios
    rad     = sqrt(-4*y[6]^2+(params["atot"]-y[6])^2);
    atp     = (params["atot"]+rad-y[6])/2;
    mgadp   = 0.165*y[6]; # magnesium adp complex
    adp3m   = 0.135*y[6]; # adp for transport into matrix  (löffler petrides s. 238)
    atp4m   = 0.05 *atp;  # atp for transport into matrix (löffler petrides s. 238)
    amp     = y[6]^2/atp; # adenosine monophosphate
    
    # flux activation functions
    topo    = 0.08+0.89*mgadp^2/params["kdd"]^2+0.16*mgadp/params["kdd"] ;
    bottomo = (1+mgadp/params["kdd"])^2*(1+atp4m/params["ktt"] + adp3m/params["ktd"]) ;
    katpo   = topo/bottomo;                                                 # I_K(ATP) activation function
    minf    = 1/(1+exp((params["vm"]-y[1])/params["sm"]));                  # ca pmca activation function
    ninf    = 1/(1+exp((params["nin"]-y[1])/params["sn"]));                 # rectifing current acti
    qinf    = y[3]^2/(params["kd"]^2+y[3]^2);                               # I_K(Ca) activation functin
    # fluxes                                                                   
    katpo = 0.009
    ik      = params["gk"]*y[2]*(y[1]-params["vk"]);                        # rectifying current
    ikca    = -params["gkca"]*qinf*(params["vk"]-y[1]);                     # ca dependent k current
    ikatp   = params["gkatpbar"]*katpo*(y[1]-params["vk"]);                 # k flux atp dependent
    ica     = params["gca"]*minf*(y[1]-params["vca"]);                      # ca flux for action potentials
    
    Jer     = params["kserca"]*y[3] - params["pleak"]*(y[4]-y[3]);          # ca flux density across er membrane
    Jm      = params["kuni"]*y[3] - params["knaca"]*(y[5]-y[3]);            # ca flux density across mitochonidium
    Jmem    = -(params["alpha"]/params["vcyt"]*ica + params["kpmca"]*y[3]); # ca flux density (cell membrane)
                                                                               
    # glycolitic oscillations                                                                                       
    Jpfk    = pfk_activity(atp, y[6], y[7], y[8], amp, params)              # pfk activity
    sinfty  = y[5]/(y[5]+params["kCaPDH"]);                                 # Michaelis-Menten function   
    Jpdh    = params["vpdh"]*sinfty*sqrt(y[8]);                             # pdh activity
    
    #save_stats(amp, atp, katpo, minf, ninf, qinf, ik, ikca, Jer, ikatp, Jm, ica, Jmem, Jpfk, sinfty, Jpdh);
    dy[1] = -(ica + ik + ikca + ikatp)/params["Cm"]
    dy[2] = -(y[2]-ninf)/params["taun"]
    dy[3] = params["fca"]*(Jmem - Jm - Jer)
    dy[4] = params["fca"]*params["sigmaer"]*Jer
    dy[5] = params["fca"]*params["sigmam"]*Jm
    dy[6] = (atp-exp((1+2.2 * Jpdh/(Jpdh+0.05)) * (1-y[3]/0.35))*y[6])/params["taua"]
    dy[7] = 0.3*(params["Jgk"]-Jpfk)
    dy[8] = Jpfk-Jpdh/2
end

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

# This function attempts to solve the ode system.
# If an DomainError exception occures it retries with
# increased accuracy. If an unexpted exception is thrown
# it is printed to stdout
# Refere to the DifferentialEquations.jl docs
# https://diffeq.sciml.ai/stable/tutorials/ode_example/
function trysolve(system, callback, iteration)
    tol = 1e-1^(iteration*6)
    backup = deepcopy(system) # Create copy for further iterations.
			      # This is necessary because the Med and ExpMed objects are destroyed 
    problem = ODEProblem(system.ode_func, system.y0, (0.0, system.time), system.params)
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
                callback = cb
        )
    catch err
        if isa(err, DomainError)
            return trysolve(backup, callback, iteration+1)
        else
            println(err)
        end
    end
end

# Takes a settings object (NamedTuple) and changes to parameters for the simulation
# accordingly. It returns DiffEq Solution object, the time series for the 8 values
# and a plot fo the solution.
function simulate(system; iteration = 1)
    if :change_params ∈ keys(system)
        for (key,val) ∈ system.change_params
                system.params[key] = val;
            end
    end
    if :meds ∈ keys(system)
        system.params["meds"] = system.meds
    end
    
    @time solution = trysolve(system, :meds ∈ keys(system), iteration)
    
    if solution == nothing
        return solution
    end
    
    matrix = scale_solution_columns(solution)
    if :plot_params ∈ keys(system) || "plot_params" ∈ keys(system)
        solution_plot = plot(
            solution.t/6000, 
            matrix[system.plot_params];
            label = hcat(labels[system.plot_params]...),
            system.plot_args...)
    else
        solution_plot = plot(solution.t/6000, matrix;
            label = hcat(labels[system.plot_params]...),
            system.plot_args...)
    end
    return (solution, matrix, solution_plot)
end

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

end # module
