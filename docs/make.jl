using PowerSeries
using Documenter

DocMeta.setdocmeta!(PowerSeries, :DocTestSetup, :(using PowerSeries); recursive = true)

makedocs(;
    modules = [PowerSeries],
    authors = "Michael Boyle <michael.oliver.boyle@gmail.com> and contributors",
    repo = "https://github.com/moble/PowerSeries.jl/blob/{commit}{path}#{line}",
    sitename = "PowerSeries.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://moble.github.io/PowerSeries.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/moble/PowerSeries.jl", devbranch = "main")
