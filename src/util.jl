## Mode function

function Φ(i::Integer, x̃::Float64)
    return sin(i * π * x̃)
end

function ∂ₓΦ(i::Integer, x̃::Float64)
    return i * π * cos(i * π * x̃)
end

function ∂ₓₓΦ(i::Integer, x̃::Float64)
    return -(i * π)^2 * sin(i * π * x̃)
end

function ∫Φ(i::Integer)
    return (1 + (-1)^(1 + i)) / (i * pi)
end

function ∫ΦΦ(i::Integer, j::Integer)
    if i == j
        return 1 / 2
    else
        return 0
    end
end

function ∫∂ₓΦ∂ₓΦ(i::Integer, j::Integer)
    if i == j
        return 1 / 2 * (i * pi)^2
    else
        return 0
    end
end

function ∫x∂ₓΦ(i::Integer)
    return (-1 + (-1)^(i)) / (i * pi)
end

function ∫x²∂ₓₓΦ(i::Integer)
    return (2 - 2 * (-1)^i + (-1)^i * i^2 * pi^2) / (i * pi)
end

function ∫x∂ₓΦΦ(i::Integer, j::Integer)
    if i == j
        return -1 / 4
    else
        return ((-1)^(i + j) * i * j) / (i^2 - j^2)
    end
end

function ∫x²∂ₓₓΦΦ(i::Integer, j::Integer)
    if i == j
        return 1 / 4 - i^2 * pi^2 / 6
    else
        return -(4 * (-1)^(i + j) * i^3 * j) / ((i^2 - j^2)^2)
    end
end


## u function

# get u(l,t) = ∑ Φᵢ(l) xᵢ
function u(x::Vector{Float64}, l::Real)
    xe = x[end]
    _u = 0

    for i = 1:N
        _u += Φ(i, l / xe) * x[i]
    end

    return _u
end

function ∂ₓu(x::Vector{Float64}, l::Real)
    xe = x[end]
    dx_u = 0

    for i = 1:N
        dx_u += ∂ₓΦ(i, l / xe) * x[i] / xe
    end

    return dx_u
end

function ∂ₜu(x::Vector{Float64}, dx::Vector{Float64}, l::Real)
    _∂ₜu = 0
    xe = x[end]
    dxe = dx[end]
    x̃ = 1 - l / xe

    for i in 1:N
        _∂ₜu += Φ(i, x̃) * dx[i] - x̃ * dxe / xe * ∂ₓΦ(i, x̃) * x[i]
    end

    return _∂ₜu
end

function ∂ₓₓu(x::Vector{Float64}, l::Real)
    xe = x[end]
    ddx_u = 0
    for i in 1:N
        ddx_u += ∂ₓₓΦ(i, l / xe) * x[i] / xe^2
    end

    return ddx_u
end

function ∂ₓₜu(x::Vector{Float64}, dx::Vector{Float64}, l::Real)
    _∂ₓₜu = 0
    xe = x[end]
    dxe = dx[end]
    x̃ = 1 - l / xe

    for i in 1:N
        _∂ₓₜu += xe * ∂ₓΦ(i, x̃) * dx[i]
        _∂ₓₜu += -dxe * ∂ₓΦ(i, x̃) * x[i]
        _∂ₓₜu += -x̃ * dxe * ∂ₓₓΦ(i, x̃) * x[i]
    end

    _∂ₓₜu = _∂ₓₜu / xe^2

    return _∂ₓₜu
end

function ∂ₓₑuₚ(pulley_id::Integer, x::Vector{Float64})
    _∂ₓₑuₚ = 0
    xe = x[end]
    for i in 1:N
        _∂ₓₑuₚ += x[i] * ∂ₓΦ(i, 1 - lp[pulley_id] / xe)
    end
    _∂ₓₑuₚ *= (lp[pulley_id] / xe^2)
    return _∂ₓₑuₚ
end

## θ function

function θ(xe::Float64)
    _θ = (-xe + L) / rd
    return _θ
end

function θp(pulley_id::Integer, x::Vector{Float64})
    xe = x[end]
    _θp = (-xe + u(x, xe - lp[pulley_id]) + L) / rp[pulley_id]
    return _θp
end

function dθ(dxe::Float64)
    _dθ = (-dxe) / rd
    return _dθ
end

function dθp(pulley_id::Integer, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]
    _∂ₜθ = (-dxe + ∂ₜu(x, dx, xe - lp[pulley_id]) + ∂ₓu(x, xe - lp[pulley_id]) * dxe) / rp[pulley_id]
    return _∂ₜθ
end

function s(x::Vector{Float64})
    xe = x[end]
    _θ = θ(xe)
    _θₚ = Vector{Float64}(undef, pulley_num)
    for i in 1:pulley_num
        _θₚ[i] = θp(i, x)
    end
    _s = [x; _θ; _θₚ]
    return _s
end

# construct ds = [dq; dθ; dθ₁, dθ₂,...,dθq]
function ds(x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _dθ = dθ(dxe)
    _dθₚ = Vector{Float64}(undef, pulley_num)
    for i in 1:pulley_num
        _dθₚ[i] = dθp(i, x, dx)
    end
    _ds = [dx; _dθ; _dθₚ]
    return _ds
end

## ζ function
function ζ(x::Vector{Float64}, dx::Vector{Float64}, x̃::Float64)
    _ζ = 0
    xe = x[end]
    dxe = dx[end]

    for i in 1:N
        _ζ += -2 * dxe / xe * x̃ * ∂ₓΦ(i, x̃) * dx[i]
        _ζ += 2 * dxe^2 / xe^2 * x̃ * ∂ₓΦ(i, x̃) * x[i]
        _ζ += dxe^2 / xe^2 * x̃^2 * ∂ₓₓΦ(i, x̃) * x[i]
    end
    return _ζ
end

function ∫ζ(x::Vector{Float64}, dx::Vector{Float64})
    _∫ζ = 0
    xe = x[end]
    dxe = dx[end]

    for i in 1:N
        _∫ζ += -2 * dxe / xe * ∫x∂ₓΦ(i) * dx[i]
        _∫ζ += 2 * dxe^2 / xe^2 * ∫x∂ₓΦ(i) * x[i]
        _∫ζ += dxe^2 / xe^2 * ∫x²∂ₓₓΦ(i) * x[i]
    end
    return _∫ζ
end

function ∫ζΦ(j::Integer, x::Vector{Float64}, dx::Vector{Float64})
    _∫ζΦ = 0
    xe = x[end]
    dxe = dx[end]

    for i in 1:N
        _∫ζΦ += -2 * dxe / xe * ∫x∂ₓΦΦ(i, j) * dx[i]
        _∫ζΦ += 2 * dxe^2 / xe^2 * ∫x∂ₓΦΦ(i, j) * x[i]
        _∫ζΦ += dxe^2 / xe^2 * ∫x²∂ₓₓΦΦ(i, j) * x[i]
    end
    return _∫ζΦ
end


# get friction force
function f(∂ₜθ::Real, T::Real, C::Real, tanh_ratio::Real=100)
    # _f = T * sign(∂ₜθ) + C * ∂ₜθ
    _f = T * tanh(tanh_ratio * ∂ₜθ) + C * ∂ₜθ
    return _f
end

function fm(dxe::Real, tanh_ratio::Real=100)
    # _f = -T * sign(dxe) - C * dxe
    _f = -Tm * tanh(tanh_ratio * dxe) - Cm * dxe
    return _f
end


function numerical_derivative(sol, t)
    ForwardDiff.derivative(t -> sol(t), t)
end

function get_cable_force_direct(x::Vector{Float64}, dx::Vector{Float64}, l::Real)
    xe = x[end]
    dxe = dx[end]
    cable_force = E * A * (∂ₓu(x, l))
    return cable_force
end

function get_mass_force_undirect(xe::Float64, ddxe::Float64)
    mass_force = k * (L - xe) - m * ddxe
    return mass_force
end