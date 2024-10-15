# cable params
const E = 137.0e9

cable_diameter = 5e-3
const A = π * (cable_diameter / 2)^2
const ρ = 1.225e-2
const L = 5.0
const a = sqrt(A * E / ρ)

# winch params
const rd = 5.0e-2
const ld = 0.0  # residual length to xe
const Id = 1e-8

# pulley params
const Ip = [1e-12]
const lp = [0.2]    # residual length to xe
const rp = [2.0e-2]
const pulley_num = length(Ip)

# mass params
const m = 1.0
const k = 100.0

# friction params
const Td = 0.0
const Cd = 0.0
const Tp = [0.0]
const Cp = [0.0]
const Tm = 0.0
const Cm = 0.0

# driven force
function T()
    return 1.0
end