module GeneralHamiltonians

# auxiliary
export Coefficient, value, isconstant
export AbelianMatrix, tompotensor

# operators, easier interface for building quantum operators incrementally, and used for TEBD. Should it really be here in this package?
export QTerm, QuantumOperator, simplify, qterms, coeff, positions, op, absorb_one_bodies, todict, apply!, shift
export InfiniteQuantumOperator, changeunitcell, oneperiod

# fermioninterface
export FermionicSymmetry, SpinfulFermionicSymmetry, ChargeCharge, SpinCharge, SpinlessFermionicSymmetry, Charge
export FermionicSymmetrySector, ChargeChargeSector, SpinChargeSector, space, ChargeSector
export twobody, fourbody, creation

using TensorKit
const TK = TensorKit
using DMRG, InfiniteDMRG
using DMRG: updateright, compute_scalartype

# auxiliary
include("auxiliary/coeff.jl")
include("auxiliary/abelianmatrix.jl")


# short range Hamiltonians
include("finite/abstractterm.jl")
include("finite/qterm.jl")
include("finite/abstractoperator.jl")
include("finite/quantumoperator.jl")
include("finite/arithmetics.jl")
include("finite/expecs.jl")
include("finite/tompo.jl")

# short range infinite Hamiltonian
include("infinite/abstractdefs.jl")
include("infinite/infinitequantumoperator.jl")
include("infinite/expecs.jl")

# utilities
include("utilities/boson_siteops.jl")
include("utilities/spin_siteops.jl")

# fermioninterface
include("utilities/fermioninterface/fermioninterface.jl")


end