

creation(pos::Int, symmetry::SpinlessFermionicSymmetry) = a_dagger(pos, sp_u1)

function spinless_fermion_site_ops()
    ph = Rep[U₁](0=>1, 1=>1)
    vacuum = oneunit(ph)
    σ₊ = zeros(vacuum ⊗ ph ← Rep[U₁](1=>1) ⊗ ph)
    copy!(block(σ₊, Irrep[U₁](1)), ones(1, 1))
    σ₋ = zeros(vacuum ⊗ ph ← Rep[U₁](-1=>1) ⊗ ph)
    copy!(block(σ₋, Irrep[U₁](0)), ones(1, 1))
    σz = ones(ph ← ph)
    copy!(block(σz, Irrep[U₁](0)), -ones(1, 1))
    JW = -σz
    return Dict("+"=>σ₊, "-"=>σ₋, "z"=>σz, "JW"=>JW)
end

const sp_u1 = spinless_fermion_site_ops()

_select_sp(symmetry::Charge) = sp_u1