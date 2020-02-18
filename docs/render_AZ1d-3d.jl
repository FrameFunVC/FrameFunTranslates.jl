
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

timingsAZ = zeros(length(fs), L,length(ds))
    errorsAZ = copy(timingsAZ)
    allocAZ = copy(timingsAZ)
    timingsAZR1 = copy(timingsAZ)
    errorsAZR1 = copy(timingsAZ)
    allocAZR1 = copy(timingsAZ)
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

    timingsAZS2 = copy(timingsAZ)
    errorsAZS2 = copy(timingsAZ)
    allocAZS2 = copy(timingsAZ)
    timingsAZS3 = copy(timingsAZ)
    errorsAZS3 = copy(timingsAZ)
    allocAZS3 = copy(timingsAZ)
    timingsAZI = copy(timingsAZ)
    errorsAZI = copy(timingsAZ)
    allocAZI = copy(timingsAZ)
    timingsI = copy(timingsAZ)
    errorsI = copy(timingsAZ)
    allocI = copy(timingsAZ)
    timingsAZSI = copy(timingsAZ)
    errorsAZSI = copy(timingsAZ)
    allocAZSI = copy(timingsAZ)
    include("fill_data.jl")

# You can leave this computation out
for (d,f) in zip(1:3,fs[1:3]), (i,n) in zip(1:10,nsds[d][1:10]), (j,p) in zip(1:3,ds[1:3]) #
    if d==1
        N = n
        Pbasis = BSplinePlatform(p)
    else
        N = ntuple(k->n,Val(d))
        Pbasis = NdBSplinePlatform(ntuple(k->p,Val(d)))
    end
    P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
    @show N, p

    # if d ==1 ||
    #         (d==2 && (n <= 129)) ||
    #         (d==3 && (p==1 && n<=28) || (p==2 && n<=22) || (p==3 && n<=22))
    #     if (i<=2)
    #         F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=AZStyle(), REG=pQR_solver, verbose=false)
    #     end
    #     F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=AZStyle(), REG=pQR_solver)
    #     errorsAZ[d,i,j] = norm(F[2]*F[4]-F[3])
    #
    #     @show timingsAZ[d,:,:]
    #     @show errorsAZ[d,:,:]
    #     @show allocAZ[d,:,:]
    # end
    # if d==1 ||
    #         (d==2 && n <= 215) ||
    #         (d==3 && ((p==1 && (n<=36) || (p==2 && (n <= 28)) || (p==3 && (n <= 28)))))
    #     if (i<=2)
    #         F, timingsAZR1[d,i,j], allocAZR1[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=ReducedAZStyle(), verbose=false)
    #     end
    #     F, timingsAZR1[d,i,j], allocAZR1[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=ReducedAZStyle())
    #     errorsAZR1[d,i,j] = norm(F[2]*F[4]-F[3])
    #
    #     @show timingsAZR1[d,:,:]
    #     @show errorsAZR1[d,:,:]
    #     @show allocAZR1[d,:,:]
    # end
    #
    if d==1
        Pbasis = CDBSplinePlatform(p)
    else
        Pbasis = NdCDBSplinePlatform(ntuple(k->p,Val(d)))
    end
    P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
    #
    # if d ==1 ||
    #        (d==2 &&((p==2 || p==3) && n <= 215) || ((p==1) &&n <= 359)) ||
    #        (d==3 &&((p==1 && n^d <=2e5) || (p==2 && n<=36) || (p==3 && n<=36)))
    #
    #     if (i<=2)
    #         F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=ReducedAZStyle(), verbose=false)
    #     end
    #     F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=ReducedAZStyle())
    #     errorsAZR2[d,i,j] = norm(F[2]*F[4]-F[3])
    #
    #     @show timingsAZR2[d,:,:]
    #     @show errorsAZR2[d,:,:]
    #     @show allocAZR2[d,:,:]
    #
    # end
    #
    # if d ==1 ||
    #        (d==2 &&((p==2 || p==3) && n <= 215) || ((p==1) &&n <= 359)) ||
    #        (d==3 &&((p==1 && n^d <=2e5) || (p==2 && n<36) || (p==3 && n<36)))
    #     if (i<=2)
    #         F, timingsAZ2[d,i,j], allocAZ2[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-8,solverstyle=AZStyle(), REG=pQR_solver, verbose=false, lraoptions=LRAOptions(atol=1e-14))
    #     end
    #     F, timingsAZ2[d,i,j], allocAZ2[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-8,solverstyle=AZStyle(), REG=pQR_solver, lraoptions=LRAOptions(atol=1e-14))
    #     errorsAZ2[d,i,j] = norm(F[2]*F[4]-F[3])
    #
    #     @show timingsAZ2[d,:,:]
    #     @show errorsAZ2[d,:,:]
    #     @show allocAZ2[d,:,:]
    # end
    #
    #
    if (d==1) || (d==2) || (d==3 && (p==1||p==2||(p==3&&n<=36))) #32
        if (i<=2)
            F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=SparseAZStyle(), verbose=false, REG=SPQR_solver)
        end
        F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=SparseAZStyle(), REG=SPQR_solver)
        errorsAZS[d,i,j] = norm(F[2]*F[4]-F[3])

        # @show timingsAZS[d,:,:]
        # @show errorsAZS[d,:,:]
        # @show allocAZS[d,:,:]
    end
    if (d==1) || (d==2) || (d==3 && (p==1||p==2||(p==3&&n<=36))) #32
        if (i <= 2)
            F, timingsS[d,i,j], allocS[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)
        end
        F, timingsS[d,i,j], allocS[d,i,j], _ = @timed approximate(f, P, N;threshold=1e-12,solverstyle=DirectStyle(), directsolver=SPQR_solver)
        errorsS[d,i,j] = norm(F[2]*F[4]-F[3])

        # @show timingsS[d,:,:]
        # @show errorsS[d,:,:]
        # @show allocS[d,:,:]
    end

    # println()
    # @show timingsAZ[d,:,:]
    # @show errorsAZ[d,:,:]
    # @show allocAZ[d,:,:]
    # @show timingsAZR1[d,:,:]
    # @show errorsAZR1[d,:,:]
    # @show allocAZR1[d,:,:]
    # @show timingsAZR2[d,:,:]
    # @show errorsAZR2[d,:,:]
    # @show allocAZR2[d,:,:]
    # @show timingsAZ2[d,:,:]
    # @show errorsAZ2[d,:,:]
    # @show allocAZ2[d,:,:]
    # @show timingsAZS[d,:,:]
    # @show errorsAZS[d,:,:]
    # @show allocAZS[d,:,:]
    # @show timingsS[d,:,:]
    # @show errorsS[d,:,:]
    # @show allocS[d,:,:]
    # println()
end


@show n1 = ns1
@show n2 = ns2.^2
@show n3 = ns3.^3
for d in 1:3
    @show timingsAZ[d,:,:]
    @show errorsAZ[d,:,:]
    @show allocAZ[d,:,:]
    @show timingsAZR1[d,:,:]
    @show errorsAZR1[d,:,:]
    @show allocAZR1[d,:,:]
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

using PGFPlotsX, LaTeXStrings, Printf, DocumentPGFPlots
ns = [n1,n2,n3]

# Timings
asymptoticAZ = [
    (1e-6n1.*log.(n1))[3:end],
    (1e-7n2.^2)[1:6],
    (1e-8n3.^((3*3-2)/3))[1:8]
]
cAZ = map(x->@sprintf("N^{%1.2f}",x), [1,2,(3*3-2)/3]);cAZ[1] *="\\log(N)"

@pgf grpoptions = {legend_pos="south east",legend_cell_align="left",legend_style={font="\\tiny"},cycle_list_name="mark list",
    group_style={group_size="3 by 1"},
        width="0.5\\textwidth", height="0.4\\textwidth",
        }

PAZ = @pgf GroupPlot(grpoptions,
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log"},
    vcat([
    [Plot(Table([ns[d][d==1 ? (3:end) : (:)],timingsAZ[d,d==1 ? (3:end) : (:),p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    Plot({style="black,dashed"},Table([ns[d][d==1 ? (3:end) : (:)][1:length(asymptoticAZ[d])],asymptoticAZ[d]])),
    LegendEntry(latexstring("\\mathcal O($(cAZ[d]))"))
    ] for d in 1:3]...)...)

asymptoticAZR1 = [
    (1e-7n1.*log.(ns1))[3:end],
    (1e-7n2.^(3/2))[1:7],
    (1e-8n3.^((3(3-1))/3))[1:9]
]
cAZR1 = map(x->@sprintf("N^{%1.2f}",x), [1,3/2,((3(3-1))/3)]);cAZR1[1] *="\\log(N)"
timingsAZR1[3,10,1] = 0
PAZR1 = @pgf GroupPlot(grpoptions,
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d][d==1 ? (3:end) : (:)],timingsAZR1[d,d==1 ? (3:end) : (:),p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d][d==1 ? (3:end) : (:)][1:length(asymptoticAZR1[d])],asymptoticAZR1[d]])),
    LegendEntry(latexstring("\\mathcal O($(cAZR1[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d][d==1 ? (3:end) : (:)],timingsAZ[d,d==1 ? (3:end) : (:),p]])) for p in 1:3]...,
    ] for d in 1:3]...)...)

asymptoticAZR2 = [
    1e-6n1[3:10],
    ((5e-7n2.^(3/2)))[3:8],
    (1e-8n3.^((3(3-1))/3))[3:10]
]
timingsAZR2[2,9,1] = 0
cAZR2 = map(x->@sprintf("N^{%1.2f}",x), [1,3/2,((3(3-1))/3)])
PAZR2 = @pgf GroupPlot(grpoptions,
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d][(3:end)],timingsAZR2[d,(3:end),p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d][(3:end)][1:length(asymptoticAZR2[d])],asymptoticAZR2[d]])),
    LegendEntry(latexstring("\\mathcal O($(cAZR2[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d][(3:end)],timingsAZR1[d,(3:end),p]])) for p in 1:3]...,
    ] for d in 1:3]...)...)

PAZ2 = @pgf GroupPlot(grpoptions,
   vcat([[
   {xlabel=latexstring("N"),xmode="log",ymode="log"},
   # Plot this plot
   vcat([
   [Plot(Table([ns[d][(3:end)],timingsAZ2[d,(3:end),p]])),
   LegendEntry(latexstring("p=$p"))] for p in 1:3
   ]...)...,
   # Add asymptotic lines
   Plot({style="black,dashed"},Table([ns[d][(3:end)][1:length(asymptoticAZR2[d])],asymptoticAZR2[d]])),
   LegendEntry(latexstring("\\mathcal O($(cAZR2[d]))")),
   # Plot previous plot
   [Plot({style="black,dotted"},Table([ns[d][(3:end)],timingsAZR2[d,(3:end),p]])) for p in 1:3]...,
   ] for d in 1:3]...)...)

asymptoticAZS = [
    1e-6n1,
    1e-5n2,
    1e-5n3
]

cAZS = map(x->@sprintf("N^{%1.2f}",x), [1,1,1])
PAZS = @pgf GroupPlot(grpoptions,
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsAZS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    # # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d],asymptoticAZS[d]])),
    LegendEntry(latexstring("\\mathcal O($(cAZS[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],timingsAZR2[d,:,p]])) for p in 1:3]...,
    ] for d in 1:3]...)...)

PAZS1 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsAZS2[d,:,p]])),
    LegendEntry(latexstring("p=$p")*", AZ")] for p in 1:3
    ]...)...,
    vcat([
    [Plot(Table([ns[d],timingsAZS3[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)

imgpath = splitdir(@__FILE__())[1]
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZtimings1d-3d"), PAZ)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZ2timings1d-3d"), PAZ2)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR1timings1d-3d"), PAZR1)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR2timings1d-3d"), PAZR2)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZStimings1d-3d"), PAZS)
# DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSAStimings1d-3d"), PAZS1)

PeAZS1 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xmode="log",ymode="log",legend_cell_align="left",legend_pos="south west"},
    vcat([
    [Plot(Table([ns[d],errorsAZS2[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    vcat([
    [Plot(Table([ns[d],errorsAZS3[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)
# DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSASerrors1d-3d"), PeAZS1)
