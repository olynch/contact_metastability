module DynamAnim
using Luxor

struct DynamicalSystem{T,S}
  init::T
  params::S
  step::Function
end

struct AnimParams
  setup::Function
  draw::Function
  frames::Int64
  steps_per_frame::Int64
  width::Int64
  height::Int64
  name::String
end

function anim(ds::DynamicalSystem{T,S},ap::AnimParams,path) where {T,S}
  sim = [zero(ds.init) for _ in 1:ap.frames]
  sim[1] = ds.init
  for i in 1:(ap.frames-1)
    sim[i+1] = copy(sim[i])
    for j in 1:ap.steps_per_frame
      ds.step(sim[i+1],ds.params)
    end
  end
  setup = ap.setup(sim[1],ap.width,ap.height)
  function scene(s,framenumber)
    ap.draw(sim[framenumber],setup)
  end
  m = Movie(ap.width,ap.height,ap.name)
  animate(m,[Scene(m,scene,1:ap.frames)],creategif=true,pathname=path)
end

end
