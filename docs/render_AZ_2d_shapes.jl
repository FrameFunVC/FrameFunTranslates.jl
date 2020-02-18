
using FrameFunTranslates, DomainSets, StaticArrays, FrameFun, LowRankApprox, StaticArrays
D1 = (0.0..0.8)^2
D2 = .35disk()+SVector(.5,.5)
D3 = .5(polardomain(x->.35(2+.5cos(5x)), (-1.0..1.0)^2)\(.15disk()+SVector(.5,0.)))+SVector(.5,.5)
D4 = .5(polardomain(x->.35(2+.5abs(cos(5.5x))), (-1.0..1.0)^2)\(.15disk()+SVector(.5,0.)))+SVector(.5,.5)


using Plots
plot(D1;size=(500,100),layout=(1,4),n=300, xlim=[0,1], ylim=[0,1])
    plot!(D2;subplot=2,n=1000, xlim=[0,1], ylim=[0,1])
    plot!(D3;subplot=3,n=1000, xlim=[0,1], ylim=[0,1])
    plot!(D4;subplot=4,n=1000, xlim=[0,1], ylim=[0,1])

using PyPlot; g = PeriodicEquispacedGrid(1000,0,1)^2
fig = figure("domains",figsize=(12,3))
subplot(141)
    PyPlot.imshow(mask(subgrid(g,D1))',aspect="equal",origin="lower",extent=(0.,1.,0.,1,),cmap="gray_r")
    subplot(142)
    PyPlot.imshow(mask(subgrid(g,D2))',aspect="equal",origin="lower",extent=(0.,1.,0.,1,),cmap="gray_r")
    subplot(143)
    PyPlot.imshow(mask(subgrid(g,D3))',aspect="equal",origin="lower",extent=(0.,1.,0.,1,),cmap="gray_r")
    subplot(144)
    PyPlot.imshow(mask(subgrid(g,D4))',aspect="equal",origin="lower",extent=(0.,1.,0.,1,),cmap="gray_r")
    PyPlot.savefig("domains",;bbox_inches="tight")


L = 10
ns = round.(Int,(10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1/2))
ds = 1:3
f = (x,y)->exp(x*y)
Ds = [D1,D2,D3,D4]


timingsAZ = zeros(length(Ds), L,length(ds))
    errorsAZ = copy(timingsAZ)
    allocAZ = copy(timingsAZ)
    timingsAZR2 = copy(timingsAZ)
    errorsAZR2 = copy(timingsAZ)
    allocAZR2 = copy(timingsAZ)
    timingsAZS = copy(timingsAZ)
    errorsAZS = copy(timingsAZ)
    allocAZS = copy(timingsAZ)
    timingsS = copy(timingsAZ)
    errorsS = copy(timingsAZ)
    allocS = copy(timingsAZ)
    include("fill_data_shapes.jl")

# You can leave this computation out
for (d,D) in zip(1:4,Ds[1:4]), (i,n) in zip(1:10,ns[1:10]), (j,p) in zip(1:3,ds[1:3]) #

    N = ntuple(k->n,Val(2))
    L = N .* 2
    Pbasis = NdCDBSplinePlatform(ntuple(k->p,Val(2)))
    P = ExtensionFramePlatform(Pbasis, D)
    @show N, p
    if (((p==2 || p==3) && n <= 215) || ((p==1) &&n <= 359))
        if (i<=2)
            F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=AZStyle(), REG=pQR_solver, verbose=false, lraoptions=LRAOptions(atol=1e-14))
        end
        F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=AZStyle(), REG=pQR_solver, lraoptions=LRAOptions(atol=1e-14))
        errorsAZ[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZ[d,:,:]
        @show errorsAZ[d,:,:]
        @show allocAZ[d,:,:]
    end

    if (((p==2 || p==3) && n <= 215) || ((p==1) &&n <= 359))
        if (i<=2)
            F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;directsolver=pQR_solver, threshold=1e-12,solverstyle=ReducedAZStyle(), verbose=false, lraoptions=LRAOptions(atol=1e-14))
        end
        F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;directsolver=pQR_solver, threshold=1e-12,solverstyle=ReducedAZStyle(), lraoptions=LRAOptions(atol=1e-14))
        errorsAZR2[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZR2[d,:,:]
        @show errorsAZR2[d,:,:]
        @show allocAZR2[d,:,:]

    end

    if (i<=2)
        F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=SparseAZStyle(), verbose=false, REG=SPQR_solver)
    end
    F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=SparseAZStyle(), REG=SPQR_solver)
    errorsAZS[d,i,j] = norm(F[2]*F[4]-F[3])

    @show timingsAZS[d,:,:]
    @show errorsAZS[d,:,:]
    @show allocAZS[d,:,:]

    if (i <= 2)
        F, timingsS[d,i,j], allocS[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)
    end
    F, timingsS[d,i,j], allocS[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=DirectStyle(), directsolver=SPQR_solver)
    errorsS[d,i,j] = norm(F[2]*F[4]-F[3])

    @show timingsS[d,:,:]
    @show errorsS[d,:,:]
    @show allocS[d,:,:]


    println()
    @show timingsAZ[d,:,:]
    @show errorsAZ[d,:,:]
    @show allocAZ[d,:,:]
    @show timingsAZR2[d,:,:]
    @show errorsAZR2[d,:,:]
    @show allocAZR2[d,:,:]
    @show timingsAZS[d,:,:]
    @show errorsAZS[d,:,:]
    @show allocAZS[d,:,:]
    @show timingsS[d,:,:]
    @show errorsS[d,:,:]
    @show allocS[d,:,:]
    println()
end

for d in 1:4
    println()
    @show timingsAZ[d,:,:]
    @show errorsAZ[d,:,:]
    @show allocAZ[d,:,:]
    @show timingsAZR2[d,:,:]
    @show errorsAZR2[d,:,:]
    @show allocAZR2[d,:,:]
    @show timingsAZS[d,:,:]
    @show errorsAZS[d,:,:]
    @show allocAZS[d,:,:]
    @show timingsS[d,:,:]
    @show errorsS[d,:,:]
    @show allocS[d,:,:]
    println()
end


using PGFPlotsX, LaTeXStrings, DocumentPGFPlots

Ptime = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
        legend_pos="north west",width="0.5\\textwidth",height="0.4\\textwidth",
        xlabel=latexstring("N"),ylabel="time (s.)",
        group_style={x_descriptions_at="edge bottom",
            y_descriptions_at="edge left", horizontal_sep="1em",
            vertical_sep="1em",group_size={"4 by 4"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([ns,timingsAZ[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:4]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns,timingsAZR2[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns,timingsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:4]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns,timingsAZS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns,timingsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:4]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns,timingsS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns,timingsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:4]...)...,
    )


Perror = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
        legend_pos="north east",width="0.5\\textwidth",height="0.4\\textwidth",
        xlabel=latexstring("N"),ylabel="residual",
        group_style={x_descriptions_at="edge bottom",
            y_descriptions_at="edge left", horizontal_sep="1em",
            vertical_sep="1em",group_size={"4 by 4"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([ns,errorsAZ[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:4]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns,errorsAZR2[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns,errorsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:4]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns,errorsAZS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns,errorsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:4]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns,errorsS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns,errorsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:4]...)...,
    )

Peff = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
    legend_pos="north east",width="0.5\\textwidth",height="0.4\\textwidth",
    xlabel="time (s.)",ylabel="residual",
    group_style={x_descriptions_at="edge bottom",
        y_descriptions_at="edge left", horizontal_sep="1em",vertical_sep="1em",
        group_size={"4 by 4"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsAZ[d,:,p],errorsAZ[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:4]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsAZR2[d,:,p],errorsAZR2[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZ[d,:,p],errorsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:4]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsAZS[d,:,p],errorsAZS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZ[d,:,p],errorsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:4]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsS[d,:,p],errorsS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZS[d,:,p],errorsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:4]...)...,
    )
imgpath = splitdir(@__FILE__())[1]
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZshapestim"), Ptime)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZshapeserr"), Perror)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZshapeseff"), Peff)
