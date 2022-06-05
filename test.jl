using IntegratedOscillatorModel

s = System(
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

s_min = System(
    time          = 140 * 6000,
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

sol,mat,pl = simulate(s);

plot(p)