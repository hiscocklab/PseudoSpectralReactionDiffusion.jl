using SafeTestsets

@time @safetestset "Solver tests" include("tests/solver_tests.jl")
@time @safetestset "Interface tests" include("tests/interface_tests.jl")
@time @safetestset "Time-dependent reaction tests" include("tests/reaction_tests.jl")
