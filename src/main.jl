## init program

# include the necessary packages
begin
    using DifferentialEquations
    using GLMakie
    using LinearAlgebra
    using DelimitedFiles
    using ForwardDiff
    using BenchmarkTools, Profile
end

# init config params
begin
    include("param.jl")
end

# simulation params
begin
    const N = 20
    tspan = (0.0, 5.0)
    tol = 1e-10
end

begin
    include("util.jl")
    include("ode.jl")
    function ode_function!(dX, X, p, t)
        x = X[1:(N+1)]
        dx = X[(N+2):end]

        _s = s(x)
        _ds = ds(x, dx)
        _M = M(_s, _ds)
        _F = F(_s, _ds)
        _G = G(x, dx)
        _C = C(x, dx)

        id_M = I(N + 1)
        zero_c = zeros(N + 1, 1)

        _A = _M * [id_M; _G]
        _B = _F - _M * [zero_c; _C]

        # id_X = I(N + 1)
        # zero_X = zeros(N + 1, N + 1)

        # LHS = [id_X zero_X; zero_X _A]
        # RHS = [dx; _B]
        # dX_tmp = LHS \ RHS
        ddx = inv(_A) * _B
        dX[1:N+1] = dx
        dX[N+2:end] = ddx
        # dX_tmp = [dx; ddx]
        # copyto!(dX, dX_tmp)
        nothing
    end
end

begin
    include("bc.jl")
    X0, dX0 = load_init_conditions()
    prob = ODEProblem(ode_function!, X0, tspan)
end

# solve problem
@time sol = solve(prob, TRBDF2(autodiff=false), reltol=tol, abstol=tol);

## post process

begin
    # plt_tspan = (0.0, 1.0)
    plt_tspan = tspan
    plot_size = 10000 * Int(plt_tspan[2] - plt_tspan[1])
    include("post.jl")
end

begin
    using MAT
    file = matopen("../data/matfile2.mat", "w")
    write(file, "sol_xe_2", sol_xe)
    write(file, "sol_dxe_2", sol_dxe)
    write(file, "sol_cable_force1_2", sol_cable_force1)
    write(file, "sol_cable_force2_2", sol_cable_force2)
    write(file, "sol_mass_force_2", sol_mass_force)

    close(file)
end