module Percolation

using Distributions,Luxor

function poisson_clock(t0::Float64,t1::Float64,rate::Float64)::Vector{Float64}
  dt = Exponential(1/rate)
  t = t0 + rand(dt)
  ts = []
  while t < t1
    push!(ts,t)
    t += rand(dt)
  end
  ts
end

@enum EventType death infectl infectr

struct Event
  type::EventType
  site::Int64
  time::Float64
end

struct Perc
  events::Vector{Event}

  function Perc(es::Vector{Event})
    new(es)
  end

  function Perc(n::Int,T::Float64,λ::Float64)
    deaths = [Event.(death, i, poisson_clock(0.,T,1.)) for i in 1:n]
    infectls = [Event.(infectl, i, poisson_clock(0.,T,λ)) for i in 2:n]
    infectrs = [Event.(infectr, i, poisson_clock(0.,T,λ)) for i in 1:(n-1)]
    events = [vcat(deaths...); vcat(infectls...); vcat(infectrs...)]
    sort!(events,by=e -> e.time)
    new(events)
  end
end

function simulate(init::BitVector, t1::Float64, perc::Perc)
  state = copy(init)
  for e in perc.events
    if e.time > t1
      break
    end
    if e.type == death
      state[e.site] = 0
    elseif e.type == infectl
      state[e.site-1] |= state[e.site]
    elseif e.type == infectr
      state[e.site+1] |= state[s.site]
    end
  end
  return state
end

function reverse_perc(perc::Perc, T::Float64)
  reverse_time_events = map(perc.events) do e
    if e.type == death
      Event(death,e.site,T - e.time)
    elseif e.type == infectl
      Event(infectr,e.site - 1,T - e.time)
    elseif e.type == infectr
      Event(infectl,e.site + 1,T - e.time)
    end
  end
  Perc(reverse(reverse_time_events))
end

const AliveInterval = Tuple{Float64, Union{Float64,Nothing}}

const AliveIntervals = Vector{AliveInterval}

struct PercResult
  ints :: Vector{AliveIntervals}
  linfects :: Vector{Vector{Float64}}
  rinfects :: Vector{Vector{Float64}}
  function PercResult(init::BitVector)
    n = length(init)
    ints = [AliveIntervals(init[i] ? [(0,nothing)] : []) for i in 1:n]
    linfects = [[] for i in 1:n]
    rinfects = [[] for i in 1:n]
    new(ints,linfects,rinfects)
  end
  function PercResult(ints::Vector{AliveIntervals},linfects::Vector{Vector{Float64}},rinfects::Vector{Vector{Float64}})
    new(ints,linfects,rinfects)
  end
end

alive(pr::PercResult, site::Int64) = length(pr.ints[site]) > 0 && pr.ints[site][end][2] == nothing

function infect(pr::PercResult,site::Int64, t::Float64)
  if !alive(pr,site)
    push!(pr.ints[site],(t,nothing))
  end
end

function kill(pr::PercResult,site::Int64,t::Float64)
  if alive(pr,site)
    pr.ints[site][end] = (pr.ints[site][end][1],t)
  end
end

function reverse_int(int::AliveInterval, T::Float64)
  if int[2] == nothing
    (0,T - int[1])
  else
    (T - int[2],T - int[1])
  end
end

function reverse_pr(pr::PercResult,T::Float64)
  PercResult(
    Vector{AliveIntervals}(map(is -> reverse(reverse_int.(is,T)), pr.ints)),
    Vector{Vector{Float64}}([[[]] ; map(ri -> reverse(T .- ri), pr.rinfects)[1:end-1]]),
    Vector{Vector{Float64}}([map(li -> reverse(T .- li), pr.linfects)[2:end] ; [[]]])
  )
end

function simulate_pr(init::BitVector, t1::Float64, perc::Perc)
  n = length(init)
  # makes the logic a bit easier if the dead ones died before we started
  pr = PercResult(init)
  for e in perc.events
    if e.time > t1
      break
    end
    if !alive(pr, e.site)
      continue
    end
    if e.type == death
      kill(pr,e.site,e.time)
    elseif e.type == infectl
      infect(pr,e.site-1,e.time)
      push!(pr.linfects[e.site],e.time)
    elseif e.type == infectr
      infect(pr,e.site+1,e.time)
      push!(pr.rinfects[e.site],e.time)
    end
  end
  pr
end

struct AnimParams
  pr::PercResult
  perc::Perc
  T::Float64
  width::Float64
  height::Float64
  frames::Int64
  function AnimParams(pr::PercResult,perc::Perc,T::Float64,w::Float64,h::Float64,fr::Int64)
    new(pr,perc,T,w,h,fr)
  end
  function AnimParams(init::BitVector, T::Float64, λ::Float64, width::Float64, height::Float64, frames::Int64)
    n = length(init)
    perc = Perc(n,T,λ)
    pr = simulate_pr(init, T, perc)
    new(pr,perc,T,width,height,frames)
  end
end

function reverse_pr(ap::AnimParams)
  n = length(ap.pr.ints)
  rev_perc = reverse_perc(ap.perc, ap.T)
  rev_init = BitVector(map(i -> alive(ap.pr,i), 1:n))
  rev_pr = reverse_pr(simulate_pr(rev_init, ap.T, rev_perc),ap.T)
  AnimParams(rev_pr,ap.perc,ap.T,ap.width,ap.height,ap.frames)
end

function draw_intervals(ap::AnimParams,cur_t::Float64)
  pr = ap.pr
  n = length(pr.ints)
  perc = ap.perc
  width = ap.width
  height = ap.height
  T = ap.T
  # height scale and width scale
  ws = width / n
  hs = height / T
  transform([0.9,0,0,-0.9,0,0])
  sethue("grey")
  pt(site,t) = Point((site - ((n+1) / 2)) * ws, (t - T/2) * hs)
  for i in 1:n
    line(pt(i,0),pt(i,T),:stroke)
    for e in perc.events
      if e.type == death
        circle(pt(e.site,e.time),4,:fill)
      elseif e.type == infectl
        arrow(pt(e.site,e.time),pt(e.site-1,e.time))
      elseif e.type == infectr
        arrow(pt(e.site,e.time),pt(e.site+1,e.time))
      end
    end
  end
  sethue("red")
  setline(4)
  for i in 1:n
    for int in pr.ints[i]
      if int[1] <= cur_t
        t2 = min(cur_t,int[2] == nothing ? T : int[2])
        line(pt(i,int[1]),pt(i,t2),:stroke)
        if t2 < min(cur_t, T)
          circle(pt(i,t2),6,:fill)
        end
      end
    end
    for t in pr.linfects[i]
      if t < cur_t
        arrow(pt(i,t),pt(i-1,t))
      end
    end
    for t in pr.rinfects[i]
      if t < cur_t
        arrow(pt(i,t),pt(i+1,t))
      end
    end
  end
end

function anim_contact_process(ap::AnimParams, path) where {T,S}
  t_scale = 1.5 * ap.T / ap.frames
  function scene(s,framenumber)
    draw_intervals(ap, t_scale * framenumber)
  end
  m = Movie(ap.width,ap.height,"Contact Process")
  animate(m,[Scene(m,scene,1:ap.frames)],creategif=true,pathname=path)
end

end
