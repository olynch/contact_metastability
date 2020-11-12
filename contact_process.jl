module CP
include("dynam_anim.jl")

using Distributions,Luxor

struct Params
  λ::Float64
end

getdefault(xs,i,def) = checkbounds(Bool,xs,i) ? @inbounds(xs[i]) : def

rate(σ::Vector{Bool},p::Params,i::Int)::Float64 = σ[i] ? 1 : p.λ * sum(getdefault(σ,i+j,0) for j in [-1,1])

rates(σ::Vector{Bool},p::Params) = [rate(σ,p,i) for i in eachindex(σ)]

function step(σ::Vector{Bool},p::Params)::Float64
  rs = rates(σ,p)
  r = sum(rs)
  if r == 0
    return Inf
  end
  p = rs/r
  i = rand(Categorical(p))
  σ[i] = !σ[i]
  rand(Exponential(1/r))
end

function setup(σ::Vector{Bool},w,h)
  Tiler(w,h,1,length(σ),margin=0)
end

function draw(σ::Vector{Bool},tiles)
  for (pos,n) in tiles
    setcolor(σ[n] ? "red" : "green")
    box(pos, tiles.tilewidth, tiles.tileheight, :fill)
  end
end

function contact_anim(init,p,frames,steps_per_frame,width,height,path)
  ds = DynamAnim.DynamicalSystem(init,p,step)
  ap = DynamAnim.AnimParams(setup,draw,frames,steps_per_frame,width,height,"Contact Process")
  DynamAnim.anim(ds,ap,path)
end

end
