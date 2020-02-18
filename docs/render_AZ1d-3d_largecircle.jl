
using FrameFunTranslates, DomainSets, StaticArrays, FrameFun, LowRankApprox

L = 10
ns1 = round.(Int,(10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1))
ns2 = round.(Int,(10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1/2))
ns3 = round.(Int,(10.0.^LinRange(log10.((1e2, 1e5))...,L)).^(1/3))
    nsds = [ns1,ns2,ns3]

ds = 1:3

f1 = (x)->exp(x)
f2 = (x,y)->exp(x*y)
f3 = (x,y,z)->exp(x*y*z)
fs = [f1,f2,f3]

timingsAZ2 = zeros(length(fs), L,length(ds))
    errorsAZ2 = copy(timingsAZ2)
    allocAZ2 = copy(timingsAZ2)
    timingsAZR2 = copy(timingsAZ2)
    errorsAZR2 = copy(timingsAZ2)
    allocAZR2 = copy(timingsAZ2)
    timingsAZS = copy(timingsAZ2)
    errorsAZS = copy(timingsAZ2)
    allocAZS = copy(timingsAZ2)
    timingsS = copy(timingsAZ2)
    errorsS = copy(timingsAZ2)
    allocS = copy(timingsAZ2)
    include("fill_data_largecircle.jl")

# You can leave this computation out
for (d,f) in zip(1:3,fs[1:3]), (i,n) in zip(1:10,nsds[d][1:10]), (j,p) in zip(1:3,ds[1:3]) #

    if d==1
        Pbasis = CDBSplinePlatform(p)
        P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
        N = n
    elseif d==2
        Pbasis = NdCDBSplinePlatform(ntuple(k->p,Val(d)))
        P = ExtensionFramePlatform(Pbasis, .4*disk() + SVector(.5,.5))
        N = (n,n)
    elseif d==3
        Pbasis = NdCDBSplinePlatform(ntuple(k->p,Val(d)))
        P = ExtensionFramePlatform(Pbasis, .4*ball() + SVector(.5,.5,.5))
        N = (n,n,n)
    end

    @show N, p
    #
    if d ==1 ||
           (d==2 &&((p==2 || p==3) && n <= 215) || ((p==1) &&n <= 359)) ||
           (d==3 &&((p==1 && n^d <=2e5) || (p==2 && n<=36) || (p==3 && n<=36)))

        if (i<=2)
            F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=ReducedAZStyle(), verbose=false)
        end
        F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=ReducedAZStyle())
        errorsAZR2[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZR2[d,:,:]
        @show errorsAZR2[d,:,:]
        @show allocAZR2[d,:,:]

    end

    if d ==1 ||
           (d==2 &&((p==2 || p==3) && n <= 215) || ((p==1) &&n <= 359)) ||
           (d==3 &&((p==1 && n^d <=2e5) || (p==2 && n<36) || (p==3 && n<36)))
        if (i<=2)
            F, timingsAZ2[d,i,j], allocAZ2[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-8,solverstyle=AZStyle(), REG=pQR_solver, verbose=false, lraoptions=LRAOptions(atol=1e-14))
        end
        F, timingsAZ2[d,i,j], allocAZ2[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-8,solverstyle=AZStyle(), REG=pQR_solver, lraoptions=LRAOptions(atol=1e-14))
        errorsAZ2[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZ2[d,:,:]
        @show errorsAZ2[d,:,:]
        @show allocAZ2[d,:,:]
    end


    if (d==1) || (d==2) || (d==3 && (p==1||p==2||(p==3&&n<=36))) #32
        if (i<=2)
            F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=SparseAZStyle(), verbose=false, REG=SPQR_solver)
        end
        F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=SparseAZStyle(), REG=SPQR_solver)
        errorsAZS[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZS[d,:,:]
        @show errorsAZS[d,:,:]
        @show allocAZS[d,:,:]
    end
    if (d==1) || (d==2) || (d==3 && (p==1||p==2||(p==3&&n<=36))) #32
        if (i <= 2)
            F, timingsS[d,i,j], allocS[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)
        end
        F, timingsS[d,i,j], allocS[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=DirectStyle(), directsolver=SPQR_solver)
        errorsS[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsS[d,:,:]
        @show errorsS[d,:,:]
        @show allocS[d,:,:]
    end

    println()
    @show timingsAZR2[d,:,:]
    @show errorsAZR2[d,:,:]
    @show allocAZR2[d,:,:]
    @show timingsAZ2[d,:,:]
    @show errorsAZ2[d,:,:]
    @show allocAZ2[d,:,:]
    @show timingsAZS[d,:,:]
    @show errorsAZS[d,:,:]
    @show allocAZS[d,:,:]
    @show timingsS[d,:,:]
    @show errorsS[d,:,:]
    @show allocS[d,:,:]
    println()
end

@show n1 = ns1
@show n2 = ns2.^2
@show n3 = ns3.^3
for d in 1:3
    @show timingsAZR2[d,:,:]
    @show errorsAZR2[d,:,:]
    @show allocAZR2[d,:,:]
    @show timingsAZ2[d,:,:]
    @show errorsAZ2[d,:,:]
    @show allocAZ2[d,:,:]
    @show timingsAZS[d,:,:]
    @show errorsAZS[d,:,:]
    @show allocAZS[d,:,:]
    @show timingsS[d,:,:]
    @show errorsS[d,:,:]
    @show allocS[d,:,:]
    println()
end



using PGFPlotsX, LaTeXStrings, DocumentPGFPlots
ns = [n1,n2,n3]
timingsAZ=timingsAZ2;errorsAZ=errorsAZ2;allocAZ=allocAZ2;
Ptime = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
        ymin=1e-3,ymax=5e2,
        legend_pos="north west",width="0.5\\textwidth",height="0.4\\textwidth",
        xlabel=latexstring("N"),ylabel="time (s.)",
        group_style={x_descriptions_at="edge bottom",
            y_descriptions_at="edge left", horizontal_sep="1em",
            vertical_sep="1em",group_size={"3 by 4"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],timingsAZ[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],timingsAZR2[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],timingsAZS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],timingsS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    )

Perror = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
        legend_pos="north east",width="0.5\\textwidth",height="0.4\\textwidth",
        ymin=1e-13,ymax=1e1,
        xlabel=latexstring("N"),ylabel="residual",
        group_style={x_descriptions_at="edge bottom",
            y_descriptions_at="edge left", horizontal_sep="1em",
            vertical_sep="1em",group_size={"3 by 4"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],errorsAZ[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],errorsAZR2[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],errorsAZS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],errorsS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    )

Peff = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
        legend_pos="north east",width="0.5\\textwidth",height="0.4\\textwidth",
        ymin=1e-13,ymax=1e1, xmin=1e-3,xmax=5e2,
        xlabel="time (s.)",ylabel="residual",
        group_style={x_descriptions_at="edge bottom",
            y_descriptions_at="edge left", horizontal_sep="1em",
            vertical_sep="1em",group_size={"3 by 4"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsAZ[d,:,p],errorsAZ[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsAZR2[d,:,p],errorsAZR2[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZ[d,:,p],errorsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsAZS[d,:,p],errorsAZS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZ[d,:,p],errorsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsS[d,:,p],errorsS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZS[d,:,p],errorsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    )

PStime = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
        ymin=1e-3,ymax=5e2,
        legend_pos="north west",width="0.5\\textwidth",height="0.4\\textwidth",
        xlabel=latexstring("N"),ylabel="time (s.)",
        group_style={x_descriptions_at="edge bottom",
            y_descriptions_at="edge left", horizontal_sep="1em",
            vertical_sep="1em",group_size={"3 by 2"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],timingsAZS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],timingsS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    )
PSerror = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
        legend_pos="north east",width="0.5\\textwidth",height="0.4\\textwidth",
        ymin=1e-14,ymax=1e9,
        xlabel=latexstring("N"),ylabel="residual",
        group_style={x_descriptions_at="edge bottom",
            y_descriptions_at="edge left", horizontal_sep="1em",
            vertical_sep="1em",group_size={"3 by 2"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],errorsAZS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],errorsS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    )
imgpath = splitdir(@__FILE__())[1]

DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSAStimings1d-3d"), PStime)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSASerrors1d-3d"), PSerror)
