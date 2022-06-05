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
                            # entry later on.)
    );
