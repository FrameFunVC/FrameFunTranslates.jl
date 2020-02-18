
using FrameFunTranslates, DomainSets, StaticArrays, FrameFun, LowRankApprox, StaticArrays

L = 10
ns2 = round.(Int,(10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1/2))
ns3 = round.(Int,(10.0.^LinRange(log10.((1e2, 1e5))...,L)).^(1/3))
    nsds = [ns2,ns3]

ds = 1:3

f2 = (x,y)->exp(x*y)
f3 = (x,y,z)->exp(x*y*z)
fs = [f2,f3]

timingsAZ = zeros(length(fs), L,length(ds))
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
    include("fill_data_circ.jl")

# You can leave this computation out
for (h,d,f) in zip(1:2,2:3,fs[1:2]), (i,n) in zip(1:10,nsds[h][1:10]), (j,p) in zip(1:3,ds[1:3]) #

    N = ntuple(k->n,Val(d))
    L = N .* 2
    Pbasis = NdCDBSplinePlatform(ntuple(k->p,Val(d)))
    if d == 2
        P = ExtensionFramePlatform(Pbasis, .35disk()+SVector(.5,.5))
    elseif d==3
        P = ExtensionFramePlatform(Pbasis, .35ball()+SVector(.5,.5,.5))
    end
    @show N, p

    if d ==1 ||
           (d==2 &&((p==2 || p==3) && n <= 215) || ((p==1) &&n <= 359)) ||
           (d==3 &&((p==1 && n^d <=2e5) || (p==2 && n<36) || (p==3 && n<36)))
        if (i<=2)
            F, timingsAZ[h,i,j], allocAZ[h,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=AZStyle(), REG=pQR_solver, verbose=false)
        end
        F, timingsAZ[h,i,j], allocAZ[h,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=AZStyle(), REG=pQR_solver)
        errorsAZ[h,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZ[h,:,:]
        @show errorsAZ[h,:,:]
        @show allocAZ[h,:,:]
    end

    if d ==1 ||
           (d==2 &&((p==2 || p==3) && n <= 215) || ((p==1) &&n <= 359)) ||
           (d==3 &&((p==1 && n^d <=2e5) || (p==2 && n<=36) || (p==3 && n<=36)))

        if (i<=2)
            F, timingsAZR2[h,i,j], allocAZR2[h,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=ReducedAZStyle(), verbose=false)
        end
        F, timingsAZR2[h,i,j], allocAZR2[h,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=ReducedAZStyle())
        errorsAZR2[h,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZR2[h,:,:]
        @show errorsAZR2[h,:,:]
        @show allocAZR2[h,:,:]

    end

    if (d==1) || (d==2) || (d==3 && (p==1||p==2||(p==3&&n<=36))) #32
        if (i<=2)
            F, timingsAZS[h,i,j], allocAZS[h,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=SparseAZStyle(), verbose=false, REG=SPQR_solver)
        end
        F, timingsAZS[h,i,j], allocAZS[h,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=SparseAZStyle(), REG=SPQR_solver)
        errorsAZS[h,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZS[h,:,:]
        @show errorsAZS[h,:,:]
        @show allocAZS[h,:,:]
    end

    if (d==1) || (d==2) || (d==3 && (p==1||p==2||(p==3&&n<=36))) #32
        if (i <= 2)
            F, timingsS[h,i,j], allocS[h,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)
        end
        F, timingsS[h,i,j], allocS[h,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=DirectStyle(), directsolver=SPQR_solver)
        errorsS[h,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsS[h,:,:]
        @show errorsS[h,:,:]
        @show allocS[h,:,:]
    end




    println()
    @show timingsAZ[h,:,:]
    @show errorsAZ[h,:,:]
    @show allocAZ[h,:,:]
    @show timingsAZR2[h,:,:]
    @show errorsAZR2[h,:,:]
    @show allocAZR2[h,:,:]
    @show timingsAZS[h,:,:]
    @show errorsAZS[h,:,:]
    @show allocAZS[h,:,:]
    @show timingsS[h,:,:]
    @show errorsS[h,:,:]
    @show allocS[h,:,:]
    println()
end

@show n2 = ns2.^2
@show n3 = ns3.^3
for h in 1:2
    @show timingsAZ[h,:,:]
    @show errorsAZ[h,:,:]
    @show allocAZ[h,:,:]
    @show timingsAZR2[h,:,:]
    @show errorsAZR2[h,:,:]
    @show allocAZR2[h,:,:]
    @show timingsAZS[h,:,:]
    @show errorsAZS[h,:,:]
    @show allocAZS[h,:,:]
    @show timingsS[h,:,:]
    @show errorsS[h,:,:]
    @show allocS[h,:,:]
    println()
end
ns = [n2,n3]
using PGFPlotsX, LaTeXStrings, DocumentPGFPlots

Ptime = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west",width="0.7\\textwidth",height="0.6\\textwidth",xlabel=latexstring("N"),ylabel="time (s.)", group_style={x_descriptions_at="edge bottom", y_descriptions_at="edge left", horizontal_sep="1em",vertical_sep="1em",group_size={"2 by 4"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[h],timingsAZ[h,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for h in 1:2]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[h],timingsAZR2[h,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[h],timingsAZ[h,:,p]])) for p in 1:3]...
    ] for h in 1:2]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[h],timingsAZS[h,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[h],timingsAZ[h,:,p]])) for p in 1:3]...
    ] for h in 1:2]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[h],timingsS[h,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[h],timingsAZS[h,:,p]])) for p in 1:3]...
    ] for h in 1:2]...)...,
    )

Perror = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
        legend_pos="north east",width="0.7\\textwidth",height="0.6\\textwidth",
        xlabel=latexstring("N"),ylabel="residual",
        group_style={x_descriptions_at="edge bottom",
            y_descriptions_at="edge left", horizontal_sep="1em",
            vertical_sep="1em",group_size={"2 by 4"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[h],errorsAZ[h,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for h in 1:2]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[h],errorsAZR2[h,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[h],errorsAZ[h,:,p]])) for p in 1:3]...
    ] for h in 1:2]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[h],errorsAZS[h,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[h],errorsAZ[h,:,p]])) for p in 1:3]...
    ] for h in 1:2]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[h],errorsS[h,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[h],errorsAZS[h,:,p]])) for p in 1:3]...
    ] for h in 1:2]...)...,
    )

Peff = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
    legend_pos="north east",width="0.7\\textwidth",height="0.6\\textwidth",
    xlabel="time (s.)",ylabel="residual",
    group_style={x_descriptions_at="edge bottom",
        y_descriptions_at="edge left", horizontal_sep="1em",vertical_sep="1em",
        group_size={"2 by 4"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsAZ[h,:,p],errorsAZ[h,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for h in 1:2]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsAZR2[h,:,p],errorsAZR2[h,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZ[h,:,p],errorsAZ[h,:,p]])) for p in 1:3]...
    ] for h in 1:2]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsAZS[h,:,p],errorsAZS[h,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZ[h,:,p],errorsAZ[h,:,p]])) for p in 1:3]...
    ] for h in 1:2]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsS[h,:,p],errorsS[h,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZS[h,:,p],errorsAZS[h,:,p]])) for p in 1:3]...
    ] for h in 1:2]...)...,
    )
