using Documenter
using BlockHaloArrays

makedocs(
    sitename = "BlockHaloArrays",
    format = Documenter.HTML(),
    modules = [BlockHaloArrays]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
