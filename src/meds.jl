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

# struct to encompasse a generic type of medication
# with potentially non linear active substance impact
@with_kw mutable struct Activa
    time :: Real
    dose :: Float64
    name :: String = ""
    fade :: Float64 = 1
    func :: Function = id
    param :: String
end

Meds = Union{Med,ExpMed,Activa}

gen_default_activa(time, dose, fade = 1.0) =
Activa(
    time = time,
    dose = dose,
    fade = fade,
    func = id,
    param = ""
)

function Dz(time, dose, fade = 1.0)
    activa = gen_default_activa(time, dose, fade)
    activa.name = "Diaxozide"
    activa.func = (x,_) -> params["gkatpbar"] + 25 * x
    activa.param = "gkatpbar"
    return activa
end

function Tolb(time, dose, fade = 1.0)
    activa = gen_default_activa(time, dose, fade)
    activa.name = "Tolbutamid"
    activa.func = (x,_) -> params["gkatpbar"] - 25 * x
    activa.param = "gkatpbar"
    return activa
end

function Tg(time, dose, fade = 1.0)
    activa = gen_default_activa(time, dose, fade)
    activa.name =  "Thapsigargin"
    activa.func = (x,_) -> params["kserca"] - 0.04x
    activa.param = "kserca"
    return activa
end

function KCl(time, dose, fade = 1.0)
    activa = gen_default_activa(time, dose, fade)
    activa.name = "KCl"
    activa.func = (x,_) -> vkf((x > 4.8) ? x : params["vk"], 130)
    activa.param = "vk"
    return activa
end

function Glucose(time, dose, fade = 1.0)
    activa = gen_default_activa(time, dose, fade)
    activa.name = "Glucose"
    activa.func = (x,_) -> x / 1000
    activa.param = "Jgk"
    return activa
end



# calculates the effects of medications
# by calling the func in the supplied activa
function handle_medication(med :: Activa, i :: Any)
    if 0 < med.dose
        i.p[med.param] = med.func(med.dose, i)
    else 
        i.p[med.param] = params[med.param]
    end
end


# legacy medication handling
function handle_medication(med :: Union{Med, ExpMed}, i::Any) 
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
            i.p["vk"] = vkf((med.dose > 4.8) ? med.dose : params["vk"], 130)
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
function calc_dosage(m::Union{Activa,ExpMed}, i)
    m.dose *= m.fade^((i.t - m.time*6000)/6000)
    handle_medication(m, i)
    return m.dose ≤ 0.01 ? true : false
end

function parse_medications(i)
    for m ∈ i.p["meds"]
        if 6000 * m.time ≤ i.t
            if calc_dosage(m, i);
                filter!(n -> n ≢ m, i.p["meds"]) # remove "empty" Activa
            end
        end
    end
end