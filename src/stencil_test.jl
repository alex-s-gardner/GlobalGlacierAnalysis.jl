using Stencils

## slope
dx(v) = ((v[7]+2*v[8]+v[9]) - (v[1]+2*v[2]+v[3]))/8;
dy(v) = ((v[3]+2*v[6]+v[9]) - (v[1]+2*v[4]+v[7]))/8;

# slope_x  = dx/lx
# slope_y  = dx/ly

## planer curvature
# https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-curvature-works.htm
ddx(v) = (v[2] + v[8]) / 2 - v[5] # / lx^2
ddy(v) = (v[4] + v[6]) / 2 - v[5] # /ly^2

#Profile_Curvature = -2(ddx + ddy) * 100 where lx/y = grid_size in x and y

# load in a raster

# crop to min bounds + some buffer

# find points in raster 

# scaling in x and y

function surface_parameters(z; lx = 1, ly = 1, ind = nothing)
    stencil = Window(1)
    A = StencilArray(z, stencil)

    if isnothing(ind)
        dzdx = mapstencil(dx, A)
        dzdy = mapstencil(dy, A)
        dzddx = mapstencil(ddx, A)
        dzddy = mapstencil(ddy, A)
    else
        dzdx = Vector{eltype(z)}(undef,length(ind));
        dzdy = copy(dzdx)
        dzddx = copy(dzdx)
        dzddy = copy(dzdx)

        for i in eachindex(ind)
            v = Stencils.stencil(A, ind[i])
            dzdx[i] = dx(v)
            dzdy[i] = dy(v)
            dzddx[i] = ddx(v)
            dzddy[i] = ddy(v)
        end
        z = z[ind]
    end

    if lx != 1
        dzdx ./= lx
        dzddx ./= lx .^ 2
    end

    if ly != 1
        dzdy ./= ly
        dzddy ./= ly.^2
    end

    return z, dzdx, dzdy, dzddx, dzddy
end;


using Stencils
n = 10000;
z = rand(n,n);

stencil = Window(1);
A = StencilArray(z, stencil);

# neighborhood functions
dx(v) = ((v[7] + 2 * v[8] + v[9]) - (v[1] + 2 * v[2] + v[3])) / 8;
dy(v) = ((v[3] + 2 * v[6] + v[9]) - (v[1] + 2 * v[4] + v[7])) / 8;
ddx(v) = (v[2] + v[8]) / 2 - v[5] # / lx^2
ddy(v) = (v[4] + v[6]) / 2 - v[5] # /ly^2

# now lets say I only need the values at k points
k = n * 200 # if k/n^2 > 0.02 then use conventional approach

ind = [(rand(1:n), rand(1:n)) for i in 1:k]
ind = CartesianIndex.(ind)

# conventional approach
@time begin
    c = similar(z)
    mapstencil!(dx, c, A);
    dzdx = c[ind]
    mapstencil!(dy, c, A)
    dzdy = c[ind]
    mapstencil!(ddx, c, A)
    dzddx = c[ind]
    mapstencil!(ddy, c, A)
    dzddy = c[ind]
end;

# point locations only
@time begin
    v = Stencils.stencil.(Ref(A), ind);
    dzdx = dx.(v)
    dzdy = dy.(v)
    dzddx = ddx.(v)
    dzddy = ddy.(v)
end;





using Stencils
using Metal

# Define a random array that the kernel will be convolved with
r = rand(10000, 10000)

# Define kernel array
k = rand(1:10,(5,5))

# Define a stencil that is the same size as the kernel array
stencil = Window(1)

# Create a stencil Kernel from the stencil and the kernel array
k = Kernel(stencil, sharpen)

# Wrap the random array and the Kernel in a StencilArray
A = StencilArray(r, k)

# use `mapstencil` with the `kernelproduct` funciton to convolve the Kernel with array
mapstencil(kernelproduct, A) 



using CUDA
r = CuArray(rand(1000, 1000))
