
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
D1 = 0.0..0.5
D2 = .5(polardomain(x->.4(2+.5cos(5x)), (-1.0..1.0)^2)\(.15disk()+SVector(.5,0.)))+SVector(.5,.5)
D3 = .4ball() + SVector(.5,.5,.5)

timingsAZ = zeros(length(fs), L,length(ds))
    errorsAZ = copy(timingsAZ)
    allocAZ = copy(timingsAZ)
    timingsAZR2 = copy(timingsAZ)
    errorsAZR2 = copy(timingsAZ)
    allocAZR2 = copy(timingsAZ)
    timingsAZ2 = copy(timingsAZ)
    errorsAZ2 = copy(timingsAZ)
    allocAZ2 = copy(timingsAZ)
    timingsAZS = copy(timingsAZ)
    errorsAZS = copy(timingsAZ)
    allocAZS = copy(timingsAZ)
    timingsS = copy(timingsAZ)
    errorsS = copy(timingsAZ)
    allocS = copy(timingsAZ)

    timingsAZI = copy(timingsAZ)
    errorsAZI = copy(timingsAZ)
    allocAZI = copy(timingsAZ)
    timingsAZSI = copy(timingsAZ)
    errorsAZSI = copy(timingsAZ)
    allocAZSI = copy(timingsAZ)
    timingsI = copy(timingsAZ)
    errorsI = copy(timingsAZ)
    allocI = copy(timingsAZ)
    timingsAZSLR = copy(timingsAZ)
    errorsAZSLR = copy(timingsAZ)
    allocAZSLR = copy(timingsAZ)
    timingsLR = copy(timingsAZ)
    errorsLR = copy(timingsAZ)
    allocLR = copy(timingsAZ)
    include("fill_lowprecision.jl")
threshold = 1e-4
for (d,f) in zip(1:3,fs[1:3]), (i,n) in zip(1:10,nsds[d][1:10]), (j,p) in zip(1:3,ds[1:3]) #
    if d==1
        N = n
        Pbasis = CDBSplinePlatform(p)
        P = ExtensionFramePlatform(Pbasis, D1)
    elseif d==2
        N = ntuple(k->n,Val(d))
        Pbasis = NdCDBSplinePlatform(ntuple(k->p,Val(d)))
        P = ExtensionFramePlatform(Pbasis, D2)
    elseif d==3
        N = ntuple(k->n,Val(d))
        Pbasis = NdCDBSplinePlatform(ntuple(k->p,Val(d)))
        P = ExtensionFramePlatform(Pbasis, D3)
    end
    @show N, p
    if (d==1) ||
        (d==2 && n<=256) ||
        (d==3 && n<=22)
        if (i <= 2)
            F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=AZStyle(), REG=pQR_solver, verbose=false, lraoptions=LRAOptions(atol=threshold,rtol=threshold))
        end
        F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=AZStyle(), REG=pQR_solver, lraoptions=LRAOptions(atol=threshold,rtol=threshold))
        errorsAZ[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZ[d,:,:]
        @show errorsAZ[d,:,:]
        @show allocAZ[d,:,:]
    end

    if (d==1) ||
        (d==2 && n<=256) ||
        (d==3 && n<=22)

        if (i <= 2)
            F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=ReducedAZStyle(), verbose=false, REG=pQR_solver, lraoptions=LRAOptions(atol=threshold,rtol=threshold))
        end
        F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=ReducedAZStyle(), REG=pQR_solver, lraoptions=LRAOptions(atol=threshold,rtol=threshold))
        errorsAZR2[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZR2[d,:,:]
        @show errorsAZR2[d,:,:]
        @show allocAZR2[d,:,:]
    end

    if (i <= 2)
        F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=SparseAZStyle(), verbose=false, REG=SPQR_solver)
    end
    F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=SparseAZStyle(), REG=SPQR_solver)
    errorsAZS[d,i,j] = norm(F[2]*F[4]-F[3])

    @show timingsAZS[d,:,:]
    @show errorsAZS[d,:,:]
    @show allocAZS[d,:,:]
    if (i <= 2)
        F, timingsS[d,i,j], allocS[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)
    end
    F, timingsS[d,i,j], allocS[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=DirectStyle(), directsolver=SPQR_solver)
    errorsS[d,i,j] = norm(F[2]*F[4]-F[3])

    @show timingsS[d,:,:]
    @show errorsS[d,:,:]
    @show allocS[d,:,:]

    # Iterative
    if (d==1 && ((p==2||p==3) || (p==1&&n<1_000_000))) ||
        (d==2 && (p==1)||(p==2&&n!=599) || (p==3&&!(n==359||n==599))) ||
        (d==3 )
        if (i <= 2)
            F, timingsAZI[d,i,j], allocAZI[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=AZStyle(), verbose=false, REG=LSQR_solver,itoptions=ItOptions(atol=threshold,btol=threshold,maxiter=n^d) )
        end
        F, timingsAZI[d,i,j], allocAZI[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=AZStyle(), REG=LSQR_solver, itoptions=ItOptions(atol=threshold,btol=threshold,maxiter=n^d))
        errorsAZI[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZI[d,:,:]
        @show errorsAZI[d,:,:]
        @show allocAZI[d,:,:]
    end
    if (i <= 2)
        F, timingsAZSI[d,i,j], allocAZSI[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=SparseAZStyle(), verbose=false, REG=LSQR_solver, itoptions=ItOptions(atol=threshold,btol=threshold,maxiter=n^d))
    end
    F, timingsAZSI[d,i,j], allocAZSI[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=SparseAZStyle(), REG=LSQR_solver, itoptions=ItOptions(atol=threshold,btol=threshold,maxiter=n^d))
    errorsAZSI[d,i,j] = norm(F[2]*F[4]-F[3])

    @show timingsAZSI[d,:,:]
    @show errorsAZSI[d,:,:]
    @show allocAZSI[d,:,:]
    if (i <= 2)
        F, timingsI[d,i,j], allocI[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=IterativeStyle(), verbose=false, iterativesolver=:lsqr, itoptions=ItOptions(atol=threshold,btol=threshold,maxiter=n^d))
    end
    F, timingsI[d,i,j], allocI[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=IterativeStyle(), iterativesolver=:lsqr, itoptions=ItOptions(atol=threshold,btol=threshold,maxiter=n^d))
    errorsI[d,i,j] = norm(F[2]*F[4]-F[3])

    @show timingsI[d,:,:]
    @show errorsI[d,:,:]
    # Low rank
    if (d==1) ||
        (d==2 && ( (p==1 && n<=356) || ((p==2||p==3) && n<=256))) ||
        (d==3 && ( (p==1 && n<=28 ) || ((p==2||p==3) && n<=22 )))
        if (i <= 2)
            F, timingsAZSLR[d,i,j], allocAZSLR[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=SparseAZStyle(), verbose=false, REG=pQR_solver, lraoptions=LRAOptions(atol=threshold,rtol=threshold))
        end
        F, timingsAZSLR[d,i,j], allocAZSLR[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=SparseAZStyle(), REG=pQR_solver, lraoptions=LRAOptions(atol=threshold,rtol=threshold))
        errorsAZSLR[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZSLR[d,:,:]
        @show errorsAZSLR[d,:,:]
        @show allocAZSLR[d,:,:]
    end
    # if (d==1 && n<=5995) ||
    #     (d==2 && n<=77) ||
    #     (d==3 )
    #     if (i <= 2)
    #         F, timingsLR[d,i,j], allocLR[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=DirectStyle(), verbose=false, directsolver=pQR_solver, lraoptions=LRAOptions(atol=threshold,rtol=threshold))
    #     end
    #     F, timingsLR[d,i,j], allocLR[d,i,j], _ = @timed approximate(f, P, N;threshold=threshold,solverstyle=DirectStyle(), directsolver=pQR_solver, lraoptions=LRAOptions(atol=threshold,rtol=threshold))
    #     errorsLR[d,i,j] = norm(F[2]*F[4]-F[3])
    #
    #     @show timingsLR[d,:,:]
    #     @show errorsLR[d,:,:]
    #     @show allocLR[d,:,:]
    # end


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
    @show timingsAZI[d,:,:]
    @show errorsAZI[d,:,:]
    @show allocAZI[d,:,:]
    @show timingsAZSI[d,:,:]
    @show errorsAZSI[d,:,:]
    @show allocAZSI[d,:,:]
    @show timingsI[d,:,:]
    @show errorsI[d,:,:]
    @show allocI[d,:,:]
    @show timingsAZSLR[d,:,:]
    @show errorsAZSLR[d,:,:]
    @show allocAZSLR[d,:,:]
    # @show timingsLR[d,:,:]
    # @show errorsLR[d,:,:]
    # @show allocLR[d,:,:]
    println()
end

for d in 1:3
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
    @show timingsAZI[d,:,:]
    @show errorsAZI[d,:,:]
    @show allocAZI[d,:,:]
    @show timingsAZSI[d,:,:]
    @show errorsAZSI[d,:,:]
    @show allocAZSI[d,:,:]
    @show timingsI[d,:,:]
    @show errorsI[d,:,:]
    @show allocI[d,:,:]
    @show timingsAZSLR[d,:,:]
    @show errorsAZSLR[d,:,:]
    @show allocAZSLR[d,:,:]
    # @show timingsLR[d,:,:]
    # @show errorsLR[d,:,:]
    # @show allocLR[d,:,:]
    println()
end


@show n1 = ns1
    @show n2 = ns2.^2
    @show n3 = ns3.^3

    using PGFPlotsX, LaTeXStrings, Printf, DocumentPGFPlots
    ns = [n1,n2,n3]
@pgf ff = d-> (d==2) ? {xmax=maximum(n2)} : {}
a
time1 = @pgf GroupPlot({cycle_list_name="mark list",legend_style={font="\\tiny"},
    xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west",
    width="0.45\\textwidth",height="0.3\\textwidth",xlabel=L"N", ylabel="time (s.)",
    ymin=5e-4,ymax=5e2,group_style={x_descriptions_at="edge bottom",
        y_descriptions_at="edge left", horizontal_sep="1em",vertical_sep="1em",group_size={"3 by 6"}}},
    vcat([[
    {ff(d)...,(d!=1 ? {legend_to_name="test"} : {})...,},
    vcat([
    [Plot(Table([ns[d],timingsAZ[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {ff(d)...,legend_to_name="test",},
    vcat([
    [Plot(Table([ns[d],timingsAZR2[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test",},
    vcat([
    [Plot(Table([ns[d],timingsAZS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test",},
    vcat([
    [Plot(Table([ns[d],timingsS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test",},
    vcat([
    [Plot(Table([ns[d],timingsAZI[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test",},
    vcat([
    [Plot(Table([ns[d],timingsI[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,

    )

time2 = @pgf GroupPlot({cycle_list_name="mark list",legend_style={font="\\tiny"},
    xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west",
    width="0.45\\textwidth",height="0.3\\textwidth",xlabel=L"N", ylabel="time (s.)",
    ymin=5e-4,ymax=5e2,group_style={x_descriptions_at="edge bottom",
        y_descriptions_at="edge left", horizontal_sep="1em",vertical_sep="1em",group_size={"3 by 3"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],timingsI[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],timingsAZSI[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsI[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],timingsAZSLR[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsI[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    # vcat([[
    # {},
    # vcat([
    # [Plot(Table([ns[d],timingsLR[d,:,p]])),
    # LegendEntry(latexstring("p=$p"))] for p in 1:3
    # ]...)...,
    # [Plot({style="black,dotted"},Table([ns[d],timingsI[d,:,p]])) for p in 1:3]...
    # ] for d in 1:3]...)...,
    )




error1 = @pgf GroupPlot({cycle_list_name="mark list",legend_style={font="\\tiny"},
    xmode="log",ymode="log",legend_cell_align="left",legend_pos="north east",
    width="0.45\\textwidth",height="0.3\\textwidth",xlabel=L"N", ylabel="residual",
    ymin=1e-13,ymax=1e4,group_style={x_descriptions_at="edge bottom",
        y_descriptions_at="edge left", horizontal_sep="1em",vertical_sep="1em",group_size={"3 by 6"}}},
    vcat([[
    {ff(d)...,(d!=1 ? {legend_to_name="test"} : {})...,},
    vcat([
    [Plot(Table([ns[d],errorsAZ[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {ff(d)...,legend_to_name="test",},
    vcat([
    [Plot(Table([ns[d],errorsAZR2[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test",},
    vcat([
    [Plot(Table([ns[d],errorsAZS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZ[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test",},
    vcat([
    [Plot(Table([ns[d],errorsS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test",},
    vcat([
    [Plot(Table([ns[d],errorsAZI[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test",},
    vcat([
    [Plot(Table([ns[d],errorsI[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZS[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,

    )

error2 = @pgf GroupPlot({cycle_list_name="mark list",legend_style={font="\\tiny"},
    xmode="log",ymode="log",legend_cell_align="left",legend_pos="north east",
    width="0.45\\textwidth",height="0.3\\textwidth",xlabel=L"N", ylabel="time (s.)",
    ymin=1e-13,ymax=1e4,group_style={x_descriptions_at="edge bottom",
        y_descriptions_at="edge left", horizontal_sep="1em",vertical_sep="1em",group_size={"3 by 3"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],errorsI[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],errorsAZSI[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsI[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d],errorsAZSLR[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsI[d,:,p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    # vcat([[
    # {},
    # vcat([
    # [Plot(Table([ns[d],errorsLR[d,:,p]])),
    # LegendEntry(latexstring("p=$p"))] for p in 1:3
    # ]...)...,
    # [Plot({style="black,dotted"},Table([ns[d],errorsI[d,:,p]])) for p in 1:3]...
    # ] for d in 1:3]...)...,
    )

imgpath = splitdir(@__FILE__())[1]
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZlowprecisiontimings1d-3d"), time1)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZlowprecisionerrors1d-3d"), error1)
