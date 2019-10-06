o = out1[1]
na = o[:thetas_expect].*NaN
names = []
for n in o[:theta0][:names]
    for i=1:length(n[2])
        push!(names, n[1])
    end
end

# t_e = (hcat([o isa Dict ? o[:thetas_expect] : na for o in out1]...),
#        hcat([o isa Dict ? o[:thetas_expect] : na for o in out2]...))

t_e = (hcat([o[:thetas_expect] for o in out1 if o isa Dict]...),
       hcat([o[:thetas_expect] for o in out2 if o isa Dict]...))

t_med = (hcat([hcat(o[:thetas_quartiles]...)[3,:] for o in out1 if o isa Dict]...),
         hcat([hcat(o[:thetas_quartiles]...)[3,:] for o in out2 if o isa Dict]...))

t_mode = (hcat([o[:thetas_mode] for o in out1 if o isa Dict]...),
          hcat([o[:thetas_mode] for o in out2 if o isa Dict]...))


i = 1
ps = []
toplot = t_e
for jj=1:size(t_e[1],1)
    push!(ps, Plots.histogram(toplot[i][jj,:], lw=0, markershape=:circle, reuse=false, label=names[jj]))
#    display(ps[end])
end
Plots.plot(ps..., layout=(size(t_e[1],1),1), reuse=false)

# Looks ok for Svalbard mar-15
