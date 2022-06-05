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