# test/types.jl

using Biodiversity
using Test

filename = "files/Allium.tsv"

cstable = Biodiversity.DiscreteCSTable(filename, scpair=true)

@test Biodiversity.ispresent(cstable, "Allium rude", "Gansu") == true
@test Biodiversity.ispresent(cstable, "Allium rude", "Beijing") == false

@test length(Biodiversity.findssfromc(cstable, "Beijing")) == 19
@test length(Biodiversity.findccfroms(cstable, "Allium wallichii")) == 6
