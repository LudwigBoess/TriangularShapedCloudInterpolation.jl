"""
    This is a port of the IDL script 'TSC.pro' by Joop Schaye from Feb 1999.
    Original: https://idlastro.gsfc.nasa.gov/ftp/pro/math/tsc.pro
"""

module TriangularShapedCloudInterpolation

    export TSCInterpolation,
           get_tsc_positions

    """
        get_tsc_positions(pos::Array{<:Real}, res_elements::Array{<:Integer})

    Returns an array with interpolation positions for a given `pos` array.
    Number of interpolation positions in each dimension is given by the `res_elements` array.
    """
    function get_tsc_positions(pos::Array{<:Real}, res_elements::Array{<:Integer})

        dim = size(pos,1)
        N   = size(pos,2)

        pos_tsc = Array{eltype(pos[1]),2}(undef, dim, N)

        @inbounds for i = 1:dim
            minx = minimum(pos[i,:])
            maxx = maximum(pos[i,:]) * (1.0+1.e-6)
            dx   = -(minx - maxx) / res_elements[i]
            pos_tsc[i,:] = ( pos[i,:] .- minx ) ./ dx
        end

        return pos_tsc
    end

    """
        get_tsc_positions(pos::Array{<:Real}, res_elements::Integer)

    Helper function if resolution elements in all dimensions are the same.
    """
    get_tsc_positions(pos::Array{<:Real}, res_elements::Integer) = get_tsc_positions(pos, [res_elements, res_elements, res_elements])

    @inline function get_distances_weights(pos::Array{<:Real}, n_steps::Integer; 
                                           wraparound::Bool, isolated::Bool)

        # Coordinates of nearest grid point (ngp)
        if wraparound
            ng = floor.(pos  .+ 0.5)
        else
            ng = floor.(pos) .+ 0.5
        end

        # Distance from sample to ngp.
        dng = ng .- pos

        if wraparound
            k2 = ng
        else
            k2 = floor.(Int64, ng .- 0.5)
        end


        # Weight of ngp.
        w2 = 0.75 .- dng .* dng

        # Point before ngp.
        k1 = k2 .- 1  # Index.
        dx  = 1.0 .- dng  # Distance to sample.
        w1 = 0.5 .* (1.5 .- dx).^2  # TSC-weight.

        # Point after ngp.
        k3 = k2 .+ 1  # Index.
        dx  = 1.0 .+ dng # Distance to sample.
        w3 = 0.5 .* (1.5 .- dx).^2  # TSC-weight.

        # Periodic boundary conditions.
        bad = findall( k2 .== 0 )
        if size(bad,1) > 0
            k1[bad] .= n_steps - 1
            if isolated
                w1[bad] .= 0.0
            end
        end

        bad = findall( k2 .== (n_steps - 1) )
        if size(bad,1) > 0
            k3[bad] .= 0
            if isolated
                w3[bad] .= 0.0
            end
        end

        if wraparound
            bad = findall( k2 .== n_steps )
            if size(bad,1) > 0
                k2[bad] .= 0
                k3[bad] .= 1
            end
        end

        return k1, k2, k3, w1, w2, w3

    end


    """
        calculate_weights!(field::Array{<:Real}, tottscweight::Array{<:Real},
                                        index::Array{<:Integer}, tscweight::Array{<:Real}, 
                                        value::Array{<:Real}, Nsamples::Integer;
                                        average::Bool=false)

    Helper function to update the weights arrays.
    """
    @inline function calculate_weights!(field::Array{<:Real}, tottscweight::Array{<:Real},
                                        index::Array{<:Integer}, tscweight::Array{<:Real}, 
                                        value::Array{<:Real}, Nsamples::Integer;
                                        average::Bool=false)

        @inbounds for j = 1:Nsamples
            field[index[j]]  += tscweight[j] * value[j]
        end
        
        if average
            @inbounds for j = 1:Nsamples
                tottscweight[index[j]] += tscweight[j]
            end
        end

        return field, tottscweight
    end


    """
        find_dim_bounds(dim::Integer)

    Helper function to get loop limits for dimension loops.
    """
    @inline function find_dim_bounds(dim::Integer)

        if dim == 1
            dimx = 3
            dimy = 1
            dimz = 1
        elseif dim == 2
            dimx = 3
            dimy = 3
            dimz = 1
        elseif dim == 3
            dimx = 3
            dimy = 3
            dimz = 3
        end

        return dimx, dimy, dimz
    end


    """
        TSCInterpolation( value::Array{<:Real}, 
                          pos::Array{<:Real},        
                          res_elements::Array{<:Integer}, 
                          average::Bool=true, 
                          wraparound::Bool=false, 
                          isolated::Bool=false    )

    Runs a TSC interpolation on the `value` array based on the provided positions.
    Returns a 3D array with interpolated values.
    """
    function TSCInterpolation(  value::Array{<:Real}, 
                                pos::Array{<:Real},        
                                res_elements::Array{<:Integer};
                                average::Bool=true, 
                                wraparound::Bool=false, 
                                isolated::Bool=false    )

        Nsamples = size(value,1)
        
        # allocate arrays for indices and weights
        kx = zeros(Int64, 3, Nsamples)
        ky = zeros(Int64, 3, Nsamples)
        kz = zeros(Int64, 3, Nsamples)

        wx = ones(3, Nsamples)
        wy = ones(3, Nsamples)
        wz = ones(3, Nsamples)

        dim = size(pos,1)

        nx = res_elements[1]

        # x direction
        kx[1,:], kx[2,:], kx[3,:], wx[1,:], wx[2,:], wx[3,:] = get_distances_weights(pos[1:1,:], res_elements[1], 
                                                                                    wraparound = wraparound, 
                                                                                    isolated   = isolated)

        # y direction
        if dim > 1
            ny = res_elements[3]
            ky[1,:], ky[2,:], ky[3,:], wy[1,:], wy[2,:], wy[3,:] = get_distances_weights(pos[2:2,:], res_elements[2], 
                                                                    wraparound = wraparound, 
                                                                    isolated   = isolated)
        else
            ny = 1
        end

        # z direction
        if dim > 2
            nz = res_elements[3]

            kz[1,:], kz[2,:], kz[3,:], wz[1,:], wz[2,:], wz[3,:] = get_distances_weights(pos[3:3,:], res_elements[3], 
                                                                    wraparound = wraparound, 
                                                                    isolated   = isolated)
        else
            nz = 1
        end

        nxny = nx*ny

        # find maximum indices for dimension loop
        dimx, dimy, dimz = find_dim_bounds(dim)

        # allocate arrays for interpolated data and weighting
        field = zeros( nx * ny * nz )
        tottscweight = zeros( nx * ny * nz )

        # run interpolation
        @inbounds @fastmath for i = 1:dimx, j = 1:dimy, k = 1:dimz

            index     = @. kx[i,:] + ky[j,:] * nx + kz[k,:] * nxny + 1
            tscweight = @. wx[i,:] * wy[j,:] * wz[k,:]

            calculate_weights!( field, tottscweight,
                                index, tscweight, value, Nsamples,
                                average=average )

        end
        
        if average
            good = findall(tottscweight .> 0.0)
            field[good] ./= tottscweight[good]
        end

        return reshape(field, (nx, ny, nz))
    end

    """
        TSCInterpolation( value::Array{<:Real}, 
                          pos::Array{<:Real},        
                          res_elements::Integer;
                          average::Bool=true, 
                          wraparound::Bool=false, 
                          isolated::Bool=false  )

    Helper function to run TSC interpolation with the same number of resolution elements in all dimensions.
    """
    function TSCInterpolation(  value::Array{<:Real}, 
                                pos::Array{<:Real},        
                                res_elements::Integer;
                                average::Bool=true, 
                                wraparound::Bool=false, 
                                isolated::Bool=false    )

        dim = size(pos,1)

        res = res_elements .* ones(Int64, dim)

        return TSCInterpolation(value, pos, res, 
                                average=average,
                                wraparound=wraparound,
                                isolated=isolated )
    end

end # module
