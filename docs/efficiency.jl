

L = 10
    ns1 = round.(Int,(10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1))
    ns2 = round.(Int,(10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1/2))
    ns3 = round.(Int,(10.0.^LinRange(log10.((1e2, 1e5))...,L)).^(1/3))
    nsds = [ns1,ns2,ns3]
    n1 = ns1
    n2 = ns2.^2
    n3 = ns3.^3
    ds = 1:3
    f1 = (x)->exp(x)
    f2 = (x,y)->exp(x*y)
    f3 = (x,y,z)->exp(x*y*z)
    fs = [f1,f2,f3]
    timingsAZR2 = zeros(length(fs), L,length(ds))
    errorsAZR2 = copy(timingsAZR2)
    allocAZR2 = copy(timingsAZR2)
    timingsAZ2 = copy(timingsAZR2)
    errorsAZ2 = copy(timingsAZR2)
    allocAZ2 = copy(timingsAZR2)
    timingsAZS = copy(timingsAZR2)
    errorsAZS = copy(timingsAZR2)
    allocAZS = copy(timingsAZR2)
    timingsS = copy(timingsAZR2)
    errorsS = copy(timingsAZR2)
    allocS = copy(timingsAZR2)

    timingsAZS = copy(timingsAZR2)
    errorsAZS = copy(timingsAZR2)
    allocAZS = copy(timingsAZR2)
    timingsAZS3 = copy(timingsAZR2)
    errorsAZS3 = copy(timingsAZR2)
    allocAZS3 = copy(timingsAZR2)
    timingsAZS2 = copy(timingsAZR2)
    errorsAZS2 = copy(timingsAZR2)
    allocAZS2 = copy(timingsAZR2)
    timingsAZI = copy(timingsAZR2)
    errorsAZI = copy(timingsAZR2)
    allocAZI = copy(timingsAZR2)
    timingsAZSI = copy(timingsAZR2)
    errorsAZSI = copy(timingsAZR2)
    allocAZSI = copy(timingsAZR2)
    timingsI = copy(timingsAZR2)
    errorsI = copy(timingsAZR2)
    allocI = copy(timingsAZR2)
    include("fill_data.jl")
    using PGFPlotsX, LaTeXStrings, Printf, DocumentPGFPlots
    ns = [n1,n2,n3]


P = @pgf GroupPlot({legend_style={font="\\tiny"},cycle_list_name="mark list",
        xmode="log",ymode="log",legend_cell_align="left",
        legend_pos="north east",width="0.5\\textwidth", height="0.4\\textwidth",
        xlabel=latexstring("N"),ymin=1e-14,ymax=1e2,ylabel="residual",
        group_style={x_descriptions_at="edge bottom", y_descriptions_at="edge left",
        horizontal_sep="1em",vertical_sep="1em",group_size={"3 by 4"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d][d==1 ? (3:end) : (:)],errorsAZ2[d,d==1 ? (3:end) : (:),p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d][d==1 ? (3:end) : (:)],errorsAZR2[d,d==1 ? (3:end) : (:),p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d][d==1 ? (3:end) : (:)],errorsAZ2[d,d==1 ? (3:end) : (:),p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d][d==1 ? (3:end) : (:)],errorsAZS[d,d==1 ? (3:end) : (:),p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d][d==1 ? (3:end) : (:)],errorsAZ2[d,d==1 ? (3:end) : (:),p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([ns[d][d==1 ? (3:end) : (:)],errorsS[d,d==1 ? (3:end) : (:),p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d][d==1 ? (3:end) : (:)],errorsAZ2[d,d==1 ? (3:end) : (:),p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    )

imgpath = splitdir(@__FILE__())[1]
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZerror1d-3d"), P)


P = @pgf GroupPlot({legend_style={font="\\tiny"},cycle_list_name="mark list",
        xmode="log",ymode="log",legend_cell_align="left",
        legend_pos="north east",width="0.5\\textwidth", height="0.4\\textwidth",
        xlabel="time (s)", ylabel="residual", ymin=1e-14,ymax=1e2,xmin=5e-4,xmax=5e2,
        group_style={x_descriptions_at="edge bottom", y_descriptions_at="edge left",
        horizontal_sep="1em",vertical_sep="1em",group_size={"3 by 4"}}},
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsAZ2[d,d==1 ? (3:end) : (:),p],errorsAZ2[d,d==1 ? (3:end) : (:),p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsAZR2[d,d==1 ? (3:end) : (:),p],errorsAZR2[d,d==1 ? (3:end) : (:),p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZ2[d,d==1 ? (3:end) : (:),p],errorsAZ2[d,d==1 ? (3:end) : (:),p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsAZS[d,d==1 ? (3:end) : (:),p],errorsAZS[d,d==1 ? (3:end) : (:),p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZ2[d,d==1 ? (3:end) : (:),p],errorsAZ2[d,d==1 ? (3:end) : (:),p]])) for p in 1:3]...
    ] for d in 1:3]...)...,
    vcat([[
    {},
    vcat([
    [Plot(Table([timingsS[d,d==1 ? (3:end) : (:),p],errorsS[d,d==1 ? (3:end) : (:),p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZS[d,d==1 ? (3:end) : (:),p],errorsAZS[d,d==1 ? (3:end) : (:),p]])) for p in 1:3]...
    ] for d in 1:3]...)...)
imgpath = splitdir(@__FILE__())[1]
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZefficiency1d-3d"), P)
