# TriangularShapedCloudInterpolation.jl

```@meta
CurrentModule = TriangularShapedCloudInterpolation
DocTestSetup = quote
    using TriangularShapedCloudInterpolation
end
```

This packages is a port of the IDL script [`TSC.pro`](https://idlastro.gsfc.nasa.gov/ftp/pro/math/tsc.pro) by Joop Schaye from Feb 1999.

# Usage

This is a small guide how to use the package

## Interpolation positions

To get an array of positions you can use for TSC interpolation you can use the helper function [`get_tsc_positions`](@ref):

```julia
pos_tsc = get_tsc_positions(pos::Array{<:Real}, res_elements::Array{<:Integer})
```

Here `pos` is an Array of positions with `pos[N_dimensions, N_entries]` and `res_elements[N_dimensions]` is the number of resolution elements you want to interpolate the data with in each dimension.
You can also just provide a single integer if you want the resolution elements in all dimensions to be the same and multiple dipatch takes care of the rest.
The returned `pos_tsc` is in row-major order since that makes working with 1-3 dimensional data in the same function easier.

## TSC interpolation

To interpolate the data (e.g. density) you need to use [`TSCInterpolation`](@ref) like so:

```julia

pos     = rand(3,1_000)
density = rand(1_000)

res_elements = [20, 20, 20]

pos_tsc = get_tsc_positions(pos, res_elements)

rho_interp = TSCInterpolation( rho, pos_tsc, res_elements, 
                                    average=true)
```

# API reference

```@index
```

```@autodocs
Modules = [TriangularShapedCloudInterpolation]
```