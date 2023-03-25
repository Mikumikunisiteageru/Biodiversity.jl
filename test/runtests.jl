# test/runtests.jl

using Biodiversity
using Test
import Aqua

Aqua.test_all(Biodiversity)
