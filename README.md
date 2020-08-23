| **Documentation**                                                 | **Build Status**                                                                                | **Licence**                                                                                |
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:| :-----------------------------------------------------------------------------------------------:|
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://LudwigBoess.github.io/TriangularShapedCloudInterpolation.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LudwigBoess.github.io/TriangularShapedCloudInterpolation.jl/dev) | [![Build Status](https://travis-ci.org/LudwigBoess/TriangularShapedCloudInterpolation.jl.svg?branch=master)](https://travis-ci.org/LudwigBoess/TriangularShapedCloudInterpolation.jl) [![codecov.io](https://codecov.io/gh/LudwigBoess/TriangularShapedCloudInterpolation.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/TriangularShapedCloudInterpolation.jl?branch=master) | [![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](LICENSE.md) |

# TriangularShapedCloudInterpolation.jl

This packages is a port of the IDL script [`TSC.pro`](https://idlastro.gsfc.nasa.gov/ftp/pro/math/tsc.pro) by Joop Schaye from Feb 1999.


## Interpolation positions

To get an array of positions you can use for TSC interpolation you can use the helper function `get_tsc_positions`:

```julia
pos_tsc = get_tsc_positions(pos::Array{<:Real}, res_elements::Array{<:Integer})
```

Here `pos` is an Array of positions with `pos[N_entries, N_dimensions]` and `res_elements[N_dimensions]` is the number of resolution elements you want to interpolate the data with in each dimension.
You can also just provide a single integer if you want the resolution elements in all dimensions to be the same and multiple dipatch takes care of the rest.

## TSC interpolation

To interpolate the data (e.g. density) you need to use `TSCInterpolation` like so:

```julia

pos     = rand(1_000,3)
density = rand(1_000)

res_elements = [20, 20, 20]

pos_tsc = get_tsc_positions(pos, res_elements)

rho_interp = TSCInterpolation( rho, pos_tsc, res_elements, 
                                    average=true)
```