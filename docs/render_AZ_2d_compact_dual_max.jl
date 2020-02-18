
using FrameFunTranslates, DomainSets, StaticArrays, FrameFun, LowRankApprox, StaticArrays
m1 = Inf
    m2 = 2
    m3 = 1
    m4 = .8
    ms = [m1,m2,m3,m4]
    L = 10
    ns = round.(Int,(10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1/2))
    ds = 1:3
    f = (x,y)->exp(x*y)
    D = .5(polardomain(x->.35(2+.5cos(5x)), (-1.0..1.0)^2)\(.15disk()+SVector(.5,0.)))+SVector(.5,.5)



timingsAZ = zeros(length(ms), L,length(ds))
    errorsAZ = copy(timingsAZ)
    allocAZ = copy(timingsAZ)
    normsAZ = copy(timingsAZ)
    timingsAZR2 = copy(timingsAZ)
    errorsAZR2 = copy(timingsAZ)
    allocAZR2 = copy(timingsAZ)
    normsAZR2 = copy(timingsAZ)
    timingsAZS = copy(timingsAZ)
    errorsAZS = copy(timingsAZ)
    allocAZS = copy(timingsAZ)
    normsAZS = copy(timingsAZ)
    timingsS = copy(timingsAZ)
    errorsS = copy(timingsAZ)
    allocS = copy(timingsAZ)
    normsS = copy(timingsAZ)
    include("fill_data_compact_dual_max.jl")

# You can leave this computation out
for (d,m) in zip(1:4,ms[1:4]), (i,n) in zip(1:10,ns[1:10]), (j,p) in zip(1:3,ds[1:3]) #

    N = ntuple(k->n,Val(2))
    L = N .* 2
    Pbasis = NdCDBSplinePlatform(ntuple(k->p,Val(2)))
    P = ExtensionFramePlatform(Pbasis, D)
    @show N, p
    if (((p==2 || p==3) && n <= 215) || ((p==1) &&n <= 359))
        if (i<=2)
            F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;compact_dual_max=m,threshold=1e-12,L=L,solverstyle=AZStyle(), REG=pQR_solver, verbose=false, lraoptions=LRAOptions(atol=1e-14))
        end
        F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;compact_dual_max=m,threshold=1e-12,L=L,solverstyle=AZStyle(), REG=pQR_solver, lraoptions=LRAOptions(atol=1e-14))
        errorsAZ[d,i,j] = norm(F[2]*F[4]-F[3])
        normsAZ[d,i,j] = norm(F[4])

        @show timingsAZ[d,:,:]
        @show errorsAZ[d,:,:]
        @show allocAZ[d,:,:]
        @show normsAZ[d,:,:]
    end

    if (((p==2 || p==3) && n <= 215) || ((p==1) &&n <= 359))
        if (i<=2)
            F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;L=L,compact_dual_max=m,directsolver=pQR_solver, threshold=1e-12,solverstyle=ReducedAZStyle(), verbose=false, lraoptions=LRAOptions(atol=1e-14))
        end
        F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;L=L,compact_dual_max=m,directsolver=pQR_solver, threshold=1e-12,solverstyle=ReducedAZStyle(), lraoptions=LRAOptions(atol=1e-14))
        errorsAZR2[d,i,j] = norm(F[2]*F[4]-F[3])
        normsAZR2[d,i,j] = norm(F[4])

        @show timingsAZR2[d,:,:]
        @show errorsAZR2[d,:,:]
        @show allocAZR2[d,:,:]
        @show normsAZR2[d,:,:]
    end

    if (i<=2)
        F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;compact_dual_max=m,threshold=1e-12,L=L,solverstyle=SparseAZStyle(), verbose=false, REG=SPQR_solver)
    end
    F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;compact_dual_max=m,threshold=1e-12,L=L,solverstyle=SparseAZStyle(), REG=SPQR_solver)
    errorsAZS[d,i,j] = norm(F[2]*F[4]-F[3])
    normsAZS[d,i,j] = norm(F[4])

    @show timingsAZS[d,:,:]
    @show errorsAZS[d,:,:]
    @show allocAZS[d,:,:]
    @show normsAZS[d,:,:]

    if (i <= 2)
        F, timingsS[d,i,j], allocS[d,i,j], _ = @timed approximate(f, P, N;compact_dual_max=m,threshold=1e-12,L=L,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)
    end
    F, timingsS[d,i,j], allocS[d,i,j], _ = @timed approximate(f, P, N;compact_dual_max=m,threshold=1e-12,L=L,solverstyle=DirectStyle(), directsolver=SPQR_solver)
    errorsS[d,i,j] = norm(F[2]*F[4]-F[3])
    normsS[d,i,j] = norm(F[4])

    @show timingsS[d,:,:]
    @show errorsS[d,:,:]
    @show allocS[d,:,:]
    @show normsS[d,:,:]

    println()
    @show timingsAZ[d,:,:]
    @show errorsAZ[d,:,:]
    @show allocAZ[d,:,:]
    @show normsAZ[d,:,:]
    @show timingsAZR2[d,:,:]
    @show errorsAZR2[d,:,:]
    @show allocAZR2[d,:,:]
    @show normsAZR2[d,:,:]
    @show timingsAZS[d,:,:]
    @show errorsAZS[d,:,:]
    @show allocAZS[d,:,:]
    @show normsAZS[d,:,:]
    @show timingsS[d,:,:]
    @show errorsS[d,:,:]
    @show allocS[d,:,:]
    @show normsS[d,:,:]
    println()
end

for d in 1:4
    println()
    @show timingsAZ[d,:,:]
    @show errorsAZ[d,:,:]
    @show allocAZ[d,:,:]
    @show normsAZ[d,:,:]
    @show timingsAZR2[d,:,:]
    @show errorsAZR2[d,:,:]
    @show allocAZR2[d,:,:]
    @show normsAZR2[d,:,:]
    @show timingsAZS[d,:,:]
    @show errorsAZS[d,:,:]
    @show allocAZS[d,:,:]
    @show normsAZS[d,:,:]
    @show timingsS[d,:,:]
    @show errorsS[d,:,:]
    @show allocS[d,:,:]
    @show normsS[d,:,:]
    println()
end

using PGFPlotsX, LaTeXStrings, DocumentPGFPlots

Ptime = @pgf GroupPlot({legend_style={font="\\tiny"},xmode="log",ymode="log",legend_cell_align="left",
        ymin=1e-3,ymax=5e2,
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

Perror = @pgf GroupPlot({legend_style={font="\\tiny"},xmode="log",ymode="log",legend_cell_align="left",
        legend_pos="north east",width="0.5\\textwidth",height="0.4\\textwidth",
        ymin=1e-13,ymax=1e1,
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

Peff = @pgf GroupPlot({cycle_list_name="mark list",legend_style={font="\\tiny"},xmode="log",ymode="log",legend_cell_align="left",
    legend_pos="north east",width="0.5\\textwidth",height="0.4\\textwidth",
    xmin=1e-3,xmax=5e2,ymin=1e-13,ymax=1e1,
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
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZcompactdualtim"), Ptime)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZcompactdualerr"), Perror)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZcompactdualeff"), Peff)
