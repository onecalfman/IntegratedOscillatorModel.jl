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