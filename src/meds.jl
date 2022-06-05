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

Meds = Union{Med,ExpMed}

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