using Documenter
using AtiyahBott

makedocs(
    sitename = "AtiyahBott",
    format = Documenter.HTML(),
    modules = [AtiyahBott]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
