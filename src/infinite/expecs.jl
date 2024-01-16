# function DMRG.expectation(psiA::InfiniteMPS, m::QTerm, psiB::InfiniteMPS, envs=environments(psiA, psiB))
# 	pos = positions(m)
# 	ops = op(m)
# 	pos_end = pos[end]
# 	util = get_trivial_leg(psiA[1])
# 	envl = leftenv(envs, pos[1])
# 	envr = rightenv(envs, pos_end)
# 	@tensor hold[-3 -2; -1] := conj(psiA[pos_end][-1, 1, 2]) * envr[3, 2] * psiB[pos_end][-3, 5, 3] * ops[end][-2, 1, 4, 5] * util[4]  
# 	for j in pos_end-1:-1:pos[1]
# 		pj = findfirst(x->x==j, pos)
# 		if isnothing(pj)
# 			hold = updateright(hold, psiA[j], pj, psiB[j])
# 		else
# 			hold = updateright(hold, psiA[j], ops[pj], psiB[j])
# 		end
# 	end
# 	@tensor r = conj(util[1]) * hold[2, 1, 3] * envl[3,2]
# 	η = leading_eigenvalue(envs)
# 	nperiod = num_period(pos[1], pos_end, unitcell_size(envs))
# 	return (r * value(coeff(m))) / (envs.trace * η^nperiod)
# end
# DMRG.expectation_canonical(m::QTerm, psi::InfiniteMPS; iscanonical::Bool=false) = iscanonical ? _expectation_canonical(m, psi) : expectation(psi, m, psi)


# function DMRG.expectation(psiA::InfiniteMPS, h::Union{QuantumOperator, InfiniteQuantumOperator}, psiB::InfiniteMPS, envs=environments(psiA, psiB))
# 	r = 0.
# 	for m in qterms(h)
# 		r += expectation(psiA, m, psiB, envs)
# 	end
# 	return r
# end
# function DMRG.expectation(h::Union{QuantumOperator, InfiniteQuantumOperator}, psi::InfiniteMPS; iscanonical::Bool=false)
# 	if iscanonical
# 		r = 0.
# 		for m in qterms(h)
# 			r += expectation_canonical(m, psi)
# 		end
# 		return r
# 	else
# 		return expectation(psi, h, psi)
# 	end
# end

# function num_period(_start::Int, _stop::Int, unitcellsize::Int)
# 	_start = InfiniteDMRG.r_start(_start, unitcellsize)
# 	_stop = InfiniteDMRG.r_stop(_start, unitcellsize)
# 	return div(_stop - _start+1, unitcellsize)
# end