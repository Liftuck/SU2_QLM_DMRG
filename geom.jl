using CairoMakie, ColorSchemes, Printf

global Nx = parse(Int,ARGS[1]) # length in x-direction (must be even!)
global Ny = parse(Int,ARGS[2]) # length in y-direction from arguments
global N = Nx*Ny
@assert Nx%2==0 "Nx must be even"

function color(i) # colors for plotting
    cmap2 = to_colormap(Reverse(:thermal))
    return Makie.interpolated_getindex(cmap2, i, (0,1))
end

function save_graph_links(name) # generates heatmap of expectation values of projector 1 - |0> <0|
    sites = []
    for y=1:Ny
        for x=1:Nx
            if x%2==0
                yoff=sqrt(3)/2
            else
                yoff=0
            end
            push!(sites,(3/2*(x-1),-(sqrt(3)*(y-1)+yoff)))
        end
    end
    p_sites = Point2f[sites...]

    right = (1.0,0.0)
    up = (-0.5,sqrt(3)/2)
    down = (-0.5,-sqrt(3)/2)

    right_line = fill(right,N)
    up_line = fill(up,N)
    down_line = fill(down,N)
    for i=1:N
        right_line[i] = right_line[i] .+ sites[i]
        up_line[i] = up_line[i] .+ sites[i]
        down_line[i] = down_line[i] .+ sites[i]
    end

    p_right = Point2f[right_line...]
    p_up = Point2f[up_line...]
    p_down = Point2f[down_line...]

    c_right = fill(color(0),N)
    c_up = fill(color(0),N)
    c_down = fill(color(0),N)

    arrow_down = BezierPath([
        MoveTo(0, 0),
        CurveTo((1, 0), (1, -sqrt(3)/2), (1.5, -sqrt(3)/2))
    ])

    arrow_up = BezierPath([
        MoveTo(0, 0),
        CurveTo((1, 0), (1, sqrt(3)/2), (1.5, sqrt(3)/2))
    ])

    arrow_next_layer = BezierPath([
        MoveTo(0, 0),
        CurveTo((0.5, 0), (0.5, -sqrt(3)/4), (0, -sqrt(3)/4)),
        LineTo(-1.5*(Nx-1), -sqrt(3)/4),
        CurveTo((-1.5*(Nx-1)-0.5, -sqrt(3)/4), (-1.5*(Nx-1)-0.5, -sqrt(3)/2), (-1.5*(Nx-1), -sqrt(3)/2))
    ])

    arrow_end = BezierPath([
        MoveTo(0, 0),
        CurveTo((0.5, 0), (0.5, -sqrt(3)/4), (0, -sqrt(3)/4)),
        LineTo(-1.5*(Nx-1), -sqrt(3)/4),
        CurveTo((-1.5*(Nx-1)-0.5, -sqrt(3)/4), (-1.5*(Nx-1)-0.5, -sqrt(3)/2), (-1.5*(Nx-1), -sqrt(3)/2)),
        LineTo(-1.5*(Nx-1)+0.5, -sqrt(3)/2),
        MoveTo(-1.5*(Nx-1)+0.5-0.15, -sqrt(3)/2+0.15),
        LineTo(-1.5*(Nx-1)+0.5,      -sqrt(3)/2),
        LineTo(-1.5*(Nx-1)+0.5-0.15, -sqrt(3)/2-0.15)
    ])

    linkwidth=50
    f = Figure()

    ax = Axis(f[1,1], autolimitaspect=1, width=3/2*Nx*linkwidth, height=sqrt(3)*(Ny+1/2)*linkwidth)
    hidedecorations!(ax)
    hidespines!(ax)
    tightlimits!(ax)

    for i=1:N
        lines!(stack([p_sites[i],p_right[i]]), color=c_right[i],linewidth=10)
        lines!(stack([p_sites[i],p_up[i]]), color=c_up[i],linewidth=10)
        lines!(stack([p_sites[i],p_down[i]]), color=c_down[i],linewidth=10)
    end

    scatter!(p_sites, color="black", markersize=19, marker=:ltriangle)

    for i=1:2:(Nx*3)
        scatter!(p_sites[i], marker = arrow_down, markersize = linkwidth, color=:transparent, strokewidth = 3, strokecolor = :red)
    end
    for i=2:2:(Nx*3)
        if i%Nx == 0
            if i==Nx*3
                scatter!(p_sites[i], marker = arrow_end, markersize = linkwidth, color=:transparent, strokewidth = 3, strokecolor = :red)
            else
                scatter!(p_sites[i], marker = arrow_next_layer, markersize = linkwidth, color=:transparent, strokewidth = 3, strokecolor = :red)
            end
        else
            scatter!(p_sites[i], marker = arrow_up, markersize = linkwidth, color=:transparent, strokewidth = 3, strokecolor = :red)
        end
    end

    #Colorbar(f[1,2], limits=(0,1), colormap=Reverse(:thermal))

    resize_to_layout!(f)

    save(name*".pdf", f, pdf_version="1.4")
end

save_graph_links(@sprintf "geom_%dx%d" Nx Ny)