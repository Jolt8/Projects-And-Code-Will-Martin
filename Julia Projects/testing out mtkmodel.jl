using Unitful, DifferentialEquations
using Plots

using ModelingToolkit
using MethodOfLines
using DomainSets

@variables t x

Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

@component heat_pipe begin
    @parameters begin
        Cp 
        P 
        k
    end
    @variables begin
        T_j(..)
    end
    @equations begin
            
    end
end 

#@connector function FluidProps()