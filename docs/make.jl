using Documenter
using TriangularShapedCloudInterpolation

makedocs(
    sitename = "TriangularShapedCloudInterpolation",
    format = Documenter.HTML(),
    modules = [TriangularShapedCloudInterpolation]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
