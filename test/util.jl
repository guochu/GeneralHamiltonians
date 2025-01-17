
function boson_site_ops_u1(d::Int)
	@assert d > 1
	ph = Rep[U₁](i-1=>1 for i in 1:d)
	vacuum = oneunit(ph)
	adag = zeros(vacuum ⊗ ph ← Rep[U₁](1=>1) ⊗ ph)
	for i in 1:d-1
		copy!(block(adag, Irrep[U₁](i)), sqrt(i) * ones(1,1))
	end
	a = zeros(vacuum ⊗ ph ← Rep[U₁](-1=>1) ⊗ ph )
	for i in 1:d-1
		copy!(block(a, Irrep[U₁](i-1)), sqrt(i) * ones(1, 1))
	end
	n = zeros(ph ← ph)
	for i in 1:d-1
		copy!(block(n, Irrep[U₁](i)), i * ones(1, 1))
	end
	return Dict("+"=>adag, "-"=>a, "n"=>n)
end

function boson_site_ops_u1_2(d::Int)
	@assert d > 1
	ph = Rep[U₁](i-1=>1 for i in 1:d)
	vacuum = oneunit(ph)
	adag = zeros(vacuum ⊗ ph ← Rep[U₁](1=>1) ⊗ ph)
	for i in 1:d-1
		copy!(block(adag, Irrep[U₁](i)), sqrt(i) * ones(1,1))
	end
	a = zeros(vacuum ⊗ ph ← Rep[U₁](-1=>1) ⊗ ph )
	for i in 1:d-1
		copy!(block(a, Irrep[U₁](i-1)), sqrt(i) * ones(1, 1))
	end
	n = zeros(ph ← ph)
	for i in 1:d-1
		copy!(block(n, Irrep[U₁](i)), i * ones(1, 1))
	end
	return Dict("+"=>adag, "-"=>a, "n"=>n)
end

function spin_site_ops_u1()
    ph = Rep[U₁](0=>1, 1=>1)
    vacuum = oneunit(ph)
    σ₊ = zeros(vacuum ⊗ ph ← Rep[U₁](1=>1) ⊗ ph)
    copy!(block(σ₊, Irrep[U₁](1)), ones(1, 1))
    σ₋ = zeros(vacuum ⊗ ph ← Rep[U₁](-1=>1) ⊗ ph)
    copy!(block(σ₋, Irrep[U₁](0)), ones(1, 1))
    σz = ones(ph ← ph)
    copy!(block(σz, Irrep[U₁](0)), -ones(1, 1))
    return Dict("+"=>σ₊, "-"=>σ₋, "z"=>σz)
end

"""
	The convention is that the creation operator on the left of the annihilation operator

By convention space_l of all the operators are vacuum
"""
function spinal_fermion_site_ops_u1_su2()
	ph = Rep[U₁×SU₂]((-0.5, 0)=>1, (0.5, 0)=>1, (0, 0.5)=>1)
	bh = Rep[U₁×SU₂]((0.5, 0.5)=>1)
	vh = oneunit(ph)
	adag = zeros(Float64, vh ⊗ ph ← bh ⊗ ph)
	copy!(block(adag, Irrep[U₁](0) ⊠ Irrep[SU₂](0.5)), ones(1,1))
	copy!(block(adag, Irrep[U₁](0.5) ⊠ Irrep[SU₂](0)), sqrt(2) * ones(1,1))

	bh = Rep[U₁×SU₂]((-0.5, 0.5)=>1)
	a = zeros(Float64, vh ⊗ ph ← bh ⊗ ph)
	copy!(block(a, Irrep[U₁](0) ⊠ Irrep[SU₂](0.5)), ones(1,1))
	copy!(block(a, Irrep[U₁](-0.5) ⊠ Irrep[SU₂](0)), -sqrt(2) * ones(1,1))


	onsite_interact = zeros(Float64, ph ← ph)
	copy!(block(onsite_interact, Irrep[U₁](0.5) ⊠ Irrep[SU₂](0)), ones(1, 1))

	JW = ones(Float64, ph ← ph)
	copy!(block(JW, Irrep[U₁](0) ⊠ Irrep[SU₂](0.5)), -ones(1, 1))

	# adagJW = TensorMap(zeros, Float64, vh ⊗ ph ← bh ⊗ ph)
	# blocks(adagJW)[Irrep[U₁](0) ⊠ Irrep[SU₂](0.5)] = ones(1,1)
	# blocks(adagJW)[Irrep[U₁](0.5) ⊠ Irrep[SU₂](0)] = -sqrt(2) * ones(1,1) 

	# hund operators
	# c↑† ⊗ c↓†
	bhr = Rep[U₁×SU₂]((1, 0)=>1)
	adagadag = ones(Float64, vh ⊗ ph ← bhr ⊗ ph)

	# c↑† ⊗ c↓, this is a spin 1 sector operator!!!
	bhr = Rep[U₁×SU₂]((0, 1)=>1)
	adaga = zeros(Float64, vh ⊗ ph ← bhr ⊗ ph)
	copy!(block(adaga, Irrep[U₁](0) ⊠ Irrep[SU₂](0.5)), ones(1, 1) * (-sqrt(3) / 2))

	n = ones(Float64, ph ← ph)
	copy!(block(n, Irrep[U₁](-0.5) ⊠ Irrep[SU₂](0)), zeros(1, 1))
	copy!(block(n, Irrep[U₁](0.5) ⊠ Irrep[SU₂](0)), 2 * ones(1, 1))

	return Dict("+"=>adag, "-"=>a, "++"=>adagadag, "+-"=>adaga, "n↑n↓"=>onsite_interact, "JW"=>JW, "n"=>n)
end

function spinal_fermion_site_ops_u1_u1()
	ph = Rep[U₁×U₁]((1, 1)=>1, (1,0)=>1, (0,1)=>1, (0,0)=>1)
	vacuum = oneunit(ph)

	# adag
	adagup = zeros(Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((1,0)=>1) ⊗ ph )
	copy!(block(adagup, Irrep[U₁](1) ⊠ Irrep[U₁](0)), ones(1,1))
	copy!(block(adagup, Irrep[U₁](1) ⊠ Irrep[U₁](1)), ones(1,1))

	adagdown = zeros(Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((0,1)=>1) ⊗ ph)
	copy!(block(adagdown, Irrep[U₁](0) ⊠ Irrep[U₁](1)), ones(1,1))
	copy!(block(adagdown, Irrep[U₁](1) ⊠ Irrep[U₁](1)), -ones(1,1))

	adag = cat(adagup, adagdown, dims=3)

	# a
	aup = zeros(Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((-1,0)=>1) ⊗ ph)
	copy!(block(aup, Irrep[U₁](0) ⊠ Irrep[U₁](0)), ones(1,1))
	copy!(block(aup, Irrep[U₁](0) ⊠ Irrep[U₁](1)), ones(1,1))

	adown = zeros(Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((0,-1)=>1) ⊗ ph)
	copy!(block(adown, Irrep[U₁](0) ⊠ Irrep[U₁](0)), ones(1,1))
	copy!(block(adown, Irrep[U₁](1) ⊠ Irrep[U₁](0)), -ones(1,1))

	a = cat(aup, - adown, dims=3)

	# hund operators
	adagadag = zeros(Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((1,1)=>1) ⊗ ph)
	copy!(block(adagadag, Irrep[U₁](1) ⊠ Irrep[U₁](1)), ones(1, 1))

	# c↑† ⊗ c↓, this is a spin 1 sector operator!!!
	up = zeros(Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((1,-1)=>1) ⊗ ph)
	copy!(block(up, Irrep[U₁](1) ⊠ Irrep[U₁](0)), ones(1,1) / (-sqrt(2)))
	middle = zeros(Float64, vacuum ⊗ ph ← vacuum ⊗ ph )
	copy!(block(middle, Irrep[U₁](1) ⊠ Irrep[U₁](0)), 0.5 * ones(1,1))
	copy!(block(middle, Irrep[U₁](0) ⊠ Irrep[U₁](1)), -0.5 * ones(1,1))
	down = zeros(Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((-1,1)=>1) ⊗ ph)
	copy!(block(down, Irrep[U₁](0) ⊠ Irrep[U₁](1)), ones(1,1) / sqrt(2))
	adaga = cat(cat(up, middle, dims=3), down, dims=3)

	onsite_interact = zeros(Float64, ph ← ph)
	copy!(block(onsite_interact, Irrep[U₁](1) ⊠ Irrep[U₁](1)), ones(1,1))

	JW = ones(Float64, ph ← ph)
	copy!(block(JW, Irrep[U₁](1) ⊠ Irrep[U₁](0)), -ones(1, 1))
	copy!(block(JW, Irrep[U₁](0) ⊠ Irrep[U₁](1)), -ones(1, 1))

	occupy = ones(Float64, ph ← ph)
	copy!(block(occupy, Irrep[U₁](0) ⊠ Irrep[U₁](0)), zeros(1, 1))
	copy!(block(occupy, Irrep[U₁](1) ⊠ Irrep[U₁](1)), 2 * ones(1, 1))
	return Dict("+"=>adag, "-"=>a, "++"=>adagadag, "+-"=>adaga, "n↑n↓"=>onsite_interact, 
		"JW"=>JW, "n"=>occupy)
end

max_error(a::Vector, b::Vector) = maximum(abs.(a - b))
