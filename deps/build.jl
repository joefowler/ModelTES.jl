println("running ModelTES build.jl: "*@__FILE__)

function installed(s::String)
  try
    a=Pkg.installed(s)
    if a == nothing
      return false
    else
      return true
    end
  catch
    return false
  end
end

if !installed("ARMA")
  println("insalling ARMA")
  Pkg.clone("https://github.com/joefowler/ARMA.jl")
end
