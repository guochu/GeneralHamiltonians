
# """
# 	assume the underlying state is canonical
# """
# function _expectation_canonical(m::QTerm, psi)
# 	isstrict(m) || throw(ArgumentError("QTerm should conserve quantum number"))
# 	isconstant(m) || throw(ArgumentError("only constant QTerm allowed"))
# 	iszero(m) && return 0.
# 	pos = positions(m)
# 	ops = op(m)
# 	pos_end = pos[end]
# 	util = get_trivial_leg(psi[1])
# 	@tensor hold[-3 -2; -1] := conj(psi[pos_end][-1, 1, 2]) * psi[pos_end][-3, 4, 2] * ops[end][-2, 1, 3, 4] * util[3] 
# 	for j in pos_end-1:-1:pos[1]
# 		pj = findfirst(x->x==j, pos)
# 		if isnothing(pj)
# 			hold = updateright(hold, psi[j], pj, psi[j])
# 		else
# 			hold = updateright(hold, psi[j], ops[pj], psi[j])
# 		end
# 	end	 
# 	# s = convert(TensorMap, psi.s[pos[1]]) 
# 	s =  psi.s[pos[1]]
# 	@tensor hnew[-1; -2] := conj(s[-1, 1]) * hold[3, 2, 1] * conj(util[2]) * s[-2, 3]
# 	return tr(hnew) * value(coeff(m))
# end

function DMRG.expectation_canonical(m::QTerm, psi::AbstractMPS)
	isconstant(m) || throw(ArgumentError("only constant QTerm allowed"))
	return expectation_canonical(convert(PartialMPO, m), psi)
end

DMRG.expectation(psiA::AbstractMPS, m::QTerm, psiB::AbstractMPS, envs=environments(psiA, psiB)) = expectation(psiA, convert(PartialMPO, m), psiB, envs)
DMRG.expectation(m::QTerm, psi::AbstractMPS, envs=environments(psi, psi)) = expectation(psi, m, psi)

function DMRG.expectation(psiA::AbstractMPS, h::QuantumOperator, psiB::AbstractMPS, envs=environments(psiA, psiB))
	(length(h) <= length(psiA)) || throw(DimensionMismatch())
	r = 0.
	for m in qterms(h)
		r += expectation(psiA, m, psiB, envs)
	end
	return r
end
function DMRG.expectation_canonical(h::QuantumOperator, psi::AbstractMPS)
	r = 0.
	for m in qterms(h)
		r += expectation_canonical(m, psi)
	end
	return r
end
DMRG.expectation(h::QuantumOperator, psi::AbstractMPS, envs=environments(psi, psi)) = expectation(psi, h, psi, envs)

