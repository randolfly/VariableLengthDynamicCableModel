begin
    # GLMakie.activate!()
    dispnew(figure) = display(GLMakie.Screen(), figure)

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

    # plot xe
    fig = Figure()
    ax_xe = Axis(fig[1, 1:2], ylabel="L-xe")
    ax_dxe = Axis(fig[2, 1:2], ylabel="dxe")
    ax_ddxe = Axis(fig[3, 1:2], ylabel="ddxe")

    lines!(ax_xe, t, L .- sol_xe)
    lines!(ax_dxe, t, sol_dxe)
    lines!(ax_ddxe, t, sol_ddxe)

    dispnew(fig)

    # function eta
    fig1 = Figure()
    ax_eta = Axis(fig1[1, 1:2], ylabel="eta")
    ax_deta = Axis(fig1[2, 1:2], ylabel="deta")

    for eta_id in 1:3
        lines!(ax_eta, t, sol_eta[:, eta_id], label="eta_$eta_id")
        lines!(ax_deta, t, sol_deta[:, eta_id], label="deta_$eta_id")
    end
    axislegend(ax_deta, "Mode Shapes", position=:rt)
    # dispnew(fig1)

    # deformation plot
    fig2 = Figure()

    ax_u1 = Axis(fig2[1, 1:2], ylabel="u(xe-rp[1])")
    ax_u2 = Axis(fig2[2, 1:2], ylabel="u(L/3)")
    lines!(ax_u1, t, sol_u1)
    lines!(ax_u2, t, sol_u2)
    # dispnew(fig2)

    # force plot
    fig3 = Figure()

    ax_mf = Axis(fig3[1, 1:2], ylabel="mass force")
    # ax_cf1 = Axis(fig3[2, 1:2], ylabel="cable force(0)")
    # ax_cf2 = Axis(fig3[3, 1:2], ylabel="cable force(xe)")
    lines!(ax_mf, t, sol_mass_force)
    # lines!(ax_cf1, t, sol_cable_force1)
    # lines!(ax_cf2, t, sol_cable_force2)

    dispnew(fig3)
end