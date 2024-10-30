begin
    # resample plot points

    t = LinRange(plt_tspan..., plot_size)
    sol_xe = zeros(plot_size)
    sol_dxe = zeros(plot_size)
    sol_ddxe = zeros(plot_size)

    sol_eta = zeros(plot_size, N)
    sol_deta = zeros(plot_size, N)
    # sol_ddeta = zeros(plot_size, N)

    sol_u1 = zeros(plot_size)
    sol_u2 = zeros(plot_size)

    sol_mass_force = zeros(plot_size)
    sol_cable_force1 = zeros(plot_size)
    sol_cable_force2 = zeros(plot_size)

    # X = [dq; q]; q=[eta(1:N); xe]
    for i in 1:plot_size
        dX = numerical_derivative(sol, t[i])
        # sol_ddeta[i, :] = dX[N+2:2*N+1]
        sol_ddxe[i] = dX[end]

        x = sol(t[i])[1:N+1]
        dx = sol(t[i])[N+2:end]

        sol_eta[i, :] = x[1:N]
        sol_xe[i] = x[end]
        sol_deta[i, :] = dx[1:N]
        sol_dxe[i] = dx[end]

        sol_u1[i] = u(x, sol_xe[i] - rp[1])
        sol_u2[i] = u(x, L / 3)

        sol_mass_force[i] = get_mass_force_undirect(sol_xe[i], sol_ddxe[i])

        sol_cable_force1[i] = get_cable_force_direct(x, dx, 0.0)
        sol_cable_force2[i] = get_cable_force_direct(x, dx, sol_xe[i])
    end

end