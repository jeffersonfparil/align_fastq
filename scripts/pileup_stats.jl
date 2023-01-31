# using Pkg; Pkg.add(["StatsBase", "ProgressMeter", "DataFrames", "Plots"])
using StatsBase, ProgressMeter, DataFrames, Plots, Plots.PlotMeasures
Plots.gr()

filename = ARGS[1]
window_size = parse(Int64, ARGS[2])
number_of_chromosomes_to_include = parse(Int64, ARGS[3])
# # filename = "/data-weedomics-1/align_fastq/test/reads/test.pileup"
# # window_size = 100
# # number_of_chromosomes_to_include = 0
# filename = "/data-weedomics-1/align_fastq/test/reads/group_A.pileup"
# window_size = 100_000
# number_of_chromosomes_to_include = 7

### Extract mean and standard deviations of the minimum, maximum, mean, and standard deviations of coverage across windows
function pileup_stats(filename::String, window_size::Int64=100_000)::DataFrames.DataFrame
    file = open(filename, "r")
    vec_chromosome = []
    vec_window = []
    vec_size = []
    vec_idx_min_MODE = []
    vec_min_MEAN = []
    vec_max_MEAN = []
    vec_ave_MEAN = []
    vec_std_MEAN = []
    vec_min_SDEV = []
    vec_max_SDEV = []
    vec_ave_SDEV = []
    vec_std_SDEV = []
    vec_idx_min = []
    vec_min = []
    vec_max = []
    vec_ave = []
    vec_std = []
    window_start_position = 0
    pb_start=0; seekend(file); pb_end=position(file); seekstart(file)
    pb=ProgressMeter.Progress(pb_end, barglyphs=ProgressMeter.BarGlyphs("[=> ]"), barlen=50, color=:yellow)
    while !eof(file)
        line = split(readline(file), "\t")
        ProgressMeter.update!(pb, position(file))
        n = length(line)
        chr = line[1]
        pos = parse(Int64, line[2])
        cov = parse.(Int64, line[4:3:n])
        α = minimum(cov)
        ω = maximum(cov)
        μ = StatsBase.mean(cov)
        σ = StatsBase.std(cov); σ = isnan(σ) ? 0.0 : σ
        idx_α = σ==0.0 ? 0 : collect(1:length(cov))[cov .== minimum(cov)] ### Set idx to the 0th pool if all pools are equally covered
        if length(vec_chromosome) == 0
            push!(vec_chromosome, chr)
            push!(vec_window, 1)
            push!(vec_size, 1)
            vec_idx_min = [join(idx_α, ",")]
            vec_min = [α]
            vec_max = [ω]
            vec_ave = [μ]
            vec_std = [σ]
            window_start_position = pos
        elseif (vec_chromosome[end] == chr) & (pos < (window_start_position + window_size))
            if sum(cov) > 0
                vec_size[end] = vec_size[end] + 1
            end
            push!(vec_idx_min, join(idx_α, ","))
            push!(vec_min, α)
            push!(vec_max, ω)
            push!(vec_ave, μ)
            push!(vec_std, σ)
        else
            push!(vec_idx_min_MODE, StatsBase.mode(vec_idx_min))
            push!(vec_min_MEAN, StatsBase.mean(vec_min))
            push!(vec_max_MEAN, StatsBase.mean(vec_max))
            push!(vec_ave_MEAN, StatsBase.mean(vec_ave))
            push!(vec_std_MEAN, StatsBase.mean(vec_std))
            push!(vec_min_SDEV, StatsBase.std(vec_min))
            push!(vec_max_SDEV, StatsBase.std(vec_max))
            push!(vec_ave_SDEV, StatsBase.std(vec_ave))
            push!(vec_std_SDEV, StatsBase.std(vec_std))
            push!(vec_chromosome, chr)
            push!(vec_window, vec_window[end] + 1)
            push!(vec_size, 1)
            vec_idx_min = [join(idx_α, ",")]
            vec_min = [α]
            vec_max = [ω]
            vec_ave = [μ]
            vec_std = [σ]
            window_start_position = pos
        end
    end
    push!(vec_idx_min_MODE, StatsBase.mode(vec_idx_min))
    push!(vec_min_MEAN, StatsBase.mean(vec_min))
    push!(vec_max_MEAN, StatsBase.mean(vec_max))
    push!(vec_ave_MEAN, StatsBase.mean(vec_ave))
    push!(vec_std_MEAN, StatsBase.mean(vec_std))
    push!(vec_min_SDEV, StatsBase.std(vec_min))
    push!(vec_max_SDEV, StatsBase.std(vec_max))
    push!(vec_ave_SDEV, StatsBase.std(vec_ave))
    push!(vec_std_SDEV, StatsBase.std(vec_std))
    close(file)
    out = DataFrames.DataFrame(chr=String.(vec_chromosome),
                               window=Int64.(vec_window),
                               n=Int64.(vec_size),
                               MODE_idx_min=String.(vec_idx_min_MODE),
                               MEAN_min=Float64.(vec_min_MEAN),
                               MEAN_max=Float64.(vec_max_MEAN),
                               MEAN_mean=Float64.(vec_ave_MEAN),
                               MEAN_sdev=Float64.(vec_std_MEAN),
                               SDEV_min=Float64.(vec_min_SDEV),
                               SDEV_max=Float64.(vec_max_SDEV),
                               SDEV_mean=Float64.(vec_ave_SDEV),
                               SDEV_sdev=Float64.(vec_std_SDEV))
    return(out)
end

function fit_bestestish(x::Vector{T}, y::Vector{T}; min_poly::Float64=0.50, max_poly::Float64=7.0, step_poly::Float64=0.50)::Tuple{Vector{T}, Vector{T}, String} where T <: Number
    # min_poly=0.50; max_poly=7.0; step_poly=0.50
    n = length(x)
    @assert(n == length(y))
    idx = sortperm(x)
    x = x[idx]
    y = y[idx]
    X = ones(n)
    vec_poly_set = collect(min_poly:step_poly:max_poly)
    m = length(vec_poly_set)
    labels = []
    RMSE = []
    for i in 1:m
        # i = 1
        χ = hcat(X, x.^vec_poly_set[i])
        for j in i:m
            # j = i
            if i != j
                χ = hcat(χ, x.^vec_poly_set[j])
            end
            β̂ = y \ χ
            ŷ = χ * β̂'
            lab = i!=j ? string(vec_poly_set[i], "-", vec_poly_set[j]) : string(vec_poly_set[i])
            rmse = sqrt(mean((y .- ŷ).^2))
            # println(string(lab, ": ", rmse))
            # @show UnicodePlots.scatterplot(x, ŷ)
            # @show UnicodePlots.scatterplot(y, ŷ)
            # println("##################################")
            push!(labels, lab)
            push!(RMSE, rmse)
        end
    end
    idx = sortperm(RMSE)
    df = DataFrames.DataFrame(id=labels[idx], RMSE=RMSE[idx])
    vec_poly = parse.(Float64, split(df.id[1], "-"))
    for p in vec_poly
        X = hcat(X, x.^p)
    end
    β̂ = y \ X
    ŷ = X * β̂'
    fit_eq = string("y ~ ", round(β̂'[1], digits=2), " + ", join(string.(round.(β̂'[2:end], digits=2), "x^", vec_poly), " + "))
    return(x, ŷ, fit_eq)
end

function plot_breadth_depth(X::DataFrames.DataFrame, number_of_chromosomes_to_include::Int64, filename::String)::String
    # @time X = pileup_stats(filename, window_size)
    ### Output plot filename
    filename_output = string(join(split(filename, ".")[1:(end-1)], "."), "-coverage_stats.png")
    # filename_output = string(join(split(filename, ".")[1:(end-1)], "."), "-coverage_stats.svg")
    ### Include only the chromosomes
    if number_of_chromosomes_to_include != 0
        vec_C_chromosomes = unique(X.chr)[1:number_of_chromosomes_to_include]
        # vec_C_chromosomes_colours = repeat(collect(Plots.palette(:blues)), Int(ceil(number_of_chromosomes_to_include/2))) # colourblind-friendly
        vec_C_chromosomes_colours = repeat(collect(Plots.palette(:Set2_3)), Int(ceil(number_of_chromosomes_to_include/2))) # colourblind-friendly
        # vec_C_chromosomes_colours = repeat(collect(Plots.palette(:Dark2_3)), Int(ceil(number_of_chromosomes_to_include/2))) # colourblind-friendly
        # vec_C_chromosomes_colours = repeat(collect(Plots.palette(:Paired_3)), Int(ceil(number_of_chromosomes_to_include/2))) # colourblind-friendly
        idx = [x ∈ vec_C_chromosomes for x in X.chr]
        df = X[idx, :]
    else
        df = X
    end

    ### Find the number of pools
    f = open(filename, "r")
    n_pools = (length(split(readline(f), "\t")) - 3) / 3
    close(f)

    ### Find the median position of each chromosome
    vec_C_chromosomes_positions = []
    for i in 1:length(vec_C_chromosomes)
        idx = df.chr .== vec_C_chromosomes[i]
        push!(vec_C_chromosomes_positions, round(median(df.window[idx])))
    end

    ### PLOT 1: Breadth of coverage (Are we more or less covering each window at similar number of sites at least once?)
    l = size(df, 1)
    vec_breadth = (df.n ./ window_size) .* 100
    mean_breadth = StatsBase.mean(vec_breadth)
    p1 = Plots.plot(df.window, vec_breadth, labels="", xlabel="Genome", ylab="Breadth of coverage\nper window (%)",
                    title="Are we covering each window at similar\nnumber of sites at least once?",
                    legend=:topright, top_margin=50px, left_margin=80px, bottom_margin=50px);
    Plots.xticks!(p1, Float64.(vec_C_chromosomes_positions), string.(vec_C_chromosomes));
    for i in 1:length(vec_C_chromosomes)
        idx = df.chr .== vec_C_chromosomes[i]
        colour = vec_C_chromosomes_colours[i]
        Plots.plot!(p1, df.window[idx], vec_breadth[idx], linecolor=colour, lab="");
    end
    Plots.plot!(p1, df.window, repeat([mean_breadth], l), linecolor=:red, labels=string("Average breadth (", round(mean_breadth, digits=2), "%; ", Int64(round(mean_breadth * window_size / 100)), " per ", window_size, "bp window)"));
    # Plots.plot!(p1, size=(800, 500))
    # savefig(p1, output_p1)

    ### PLOT 2: Depth of coverage (Among the sites covered at least once, are we coverin them around the same number of times?)
    ### Include only the windows that we're covered at least once
    df = df[df.n .> 0.0, :]
    l = size(df, 1)
    ### Plot
    mean_MEAN_mean_depth = StatsBase.mean(df.MEAN_mean)
    p2 = Plots.plot(df.window, df.MEAN_mean, labels="", xlabel="Genome", ylab="Mean depth of coverage\nper window",
                    title="Among the sites covered at least once, are we\ncovering them around the same number of times?",
                    legend=:topright, top_margin=50px, left_margin=80px, bottom_margin=25px);
    Plots.xticks!(p2, Float64.(vec_C_chromosomes_positions), string.(vec_C_chromosomes));
    for i in 1:length(vec_C_chromosomes)
        idx = df.chr .== vec_C_chromosomes[i]
        colour = vec_C_chromosomes_colours[i]
        Plots.plot!(p2, df.window[idx], df.MEAN_mean[idx], linecolor=colour, lab="");
    end
    Plots.plot!(p2, df.window, repeat([mean_MEAN_mean_depth], l), linecolor=:red, labels=string("Average mean depth (", round(mean_MEAN_mean_depth, digits=2), "X)"));
    # Plots.plot!(p2, size=(800, 500))
    # savefig(p2, output_p2)

    ### PLOT 3: Mean depth of breadth of coverage with sd bars
    p3 = Plots.plot(df.window, df.MEAN_mean, yerr=df.SDEV_mean, linecolor=nothing, labels="", xlabel="Genome", ylab="Depth of coverage\n(± standard deviation)",
                    title="How much is the variation in average depth per window?",
                    legend=:topright, top_margin=50px, left_margin=80px, bottom_margin=50px);
    Plots.xticks!(p3, Float64.(vec_C_chromosomes_positions), string.(vec_C_chromosomes));
    for i in 1:length(vec_C_chromosomes)
        idx = df.chr .== vec_C_chromosomes[i]
        colour = vec_C_chromosomes_colours[i]
        Plots.plot!(p3, df.window[idx], df.MEAN_mean[idx], yerr=df.SDEV_mean[idx], linecolor=colour, markerstrokecolor=colour, labels="");
        Plots.scatter!(p3, df.window[idx], df.MEAN_mean[idx], color=:black, markerstrokewidth=0.001, markersize=2, alpha=0.50, labels="");
    end

    ### PLOT 4: Breadth and depth interactions (Is there significant trade-off between breadth and depth of coverage?)
    x = df.n .* 100 ./ window_size
    y = df.MEAN_mean
    x̂, ŷ, fit_eq = fit_bestestish(x, y)
    point_colour = Plots.palette(:tab10)[end]
    p4 = Plots.scatter(x, y, labels="", xlab="Breadth (%)", ylab="Average depth (X)",
                    title="Is there significant trade-off\nbetween breadth and depth of coverage?",
                    legend=:topright, color=point_colour, markerstrokewidth=0.001, alpha=0.25,
                    top_margin=50px, left_margin=80px, bottom_margin=50px);
    Plots.plot!(p4, x̂, ŷ, linecolor=:red, labels=fit_eq);

    ### PLOTS 5 and 6: The most frequent least covered pool
    x0 = split.(df.MODE_idx_min, ",")
    x = []
    for i in x0
        append!(x, collect(i))
    end
    x = parse.(Int, x)
    x = x[x .> 0.0] ### remove instances where all pools are equally covered
    
    ### keep the least covered pools where the coverage is exactly zero
    y0 = split.(df.MODE_idx_min[df.MEAN_min .== 0], ",")
    y = []
    for i in y0
        append!(y, collect(i))
    end
    y = parse.(Int, y)
    y = y[y .> 0.0] ### remove instances where all pools are equally covered
    p5 = Plots.histogram(x,
                         bins=1:(n_pools+1),
                         label="",
                         xlabel="Sample",
                         ylab="Frequency",
                         title="Which sample has the least coverage?",
                         top_margin=50px, left_margin=80px, bottom_margin=50px);
    xlims!(p5, 1, n_pools+1);
    Plots.xticks!(p5, collect(1:1:n_pools) .+ 0.5, string.(Int.(collect(1:1:n_pools))));
    p6 = Plots.histogram(y,
                         bins=1:(n_pools+1),
                         label="",
                         xlabel="Sample",
                         ylab="Frequency",
                         title="Which sample has the most complete lack of coverage?",
                         top_margin=50px, left_margin=80px, bottom_margin=50px);
    xlims!(p6, 1, n_pools+1);
    Plots.xticks!(p6, collect(1:1:n_pools) .+ 0.5, string.(Int.(collect(1:1:n_pools))));
    
    ### Merge plots
    pout = Plots.plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(1_600, 1_500));
    savefig(pout, filename_output)
    return(filename_output)
end

@time X = pileup_stats(filename, window_size)
@time png_out = plot_breadth_depth(X, number_of_chromosomes_to_include, filename)
