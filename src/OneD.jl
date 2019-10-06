##############################
# 1D forward model
##############################
#=
Need to become clear what is on staggered grid and what not.

Staggered:
- qtot,

Normal grid
- btilde, h, fsl


=#

# NOTE: do not access gb.fsl, gb.btilde, gb.temp, etc here!  Use
# values passed in as function arguments

"""
    fwdm1d(gb, btilde, fsl, temp, pp, pm, pn)

The 1D forward model

Returns 1D fields:
- `h` ice thickness (on nodes, i.e. length(gb.bands)+1)
- `tau_local` local basal drag (τ_b) (on bands)
- 'tau_mean` average basal drag over several ice thicknesses (on nodes)
- 'iv` surface ice velocity (on nodes)
- `qtot` total ice flux (on nodes)
"""
function fwdm1d(gb::Bands,
                # fields defined on bands:
                btilde::AbstractVector,
                fsl::AbstractVector,
                temp::AbstractVector,
                pp::Phys, pm::MPara, pn::Num)
    qtot_n = flux_n(gb, btilde) # [m^3/s]
    if pn.error_on_neg_qtot &&  any(qtot.<0)
        ii = find(qtot.<0)
        error("qtot at band indices $ii is negative!")
    end

    hs_n, tau_local, niter = calc_thick(gb, qtot_n, fsl, temp, pp, pm, pn)
    tau_mean_n  = F[bands_mean_fn(tau_local, hs_n, gb, i_n, pm.mean_tau_dist)
                    for i_n=1:length(gb)+1]
    # iv_n = calc_iv(gb.gl, gb, fsl, temp, hs_n, tau_mean_n, pp, pm, qtot_n)
    iv_n = calc_iv_masscons(gb.gl, gb, fsl, qtot_n, hs_n, pp, pm)

    return hs_n, tau_local, tau_mean_n, iv_n, qtot_n
end

"""
$(SIGNATURES)

Calculates ice surface velocity using H&F eq 4.
"""
function calc_iv(gl::Glacier, gb::Bands, fsl, temp, hs_n, tau_mean_n, pp::Phys, pm::MPara, qtot_n)
    @unpack n, rho, g, r = pp
    @unpack bands, ele, ws = gb
    A = ice_flow_factor.(on_nodes(temp))
    fsl_n = on_nodes(fsl)

    # Depth-av flow with mass-cons
    # ws_n = on_nodes(ws)
    # Plots.plot(qtot_n ./ (ws_n.*hs_n), reuse=false)
    # versus depth-av flow Paterson 1994 eq 11.22
    # display(Plots.plot!(2A./(n+2) .* (tau_mean_n./hs_n).^n.*hs_n.^(n+1) ))
    # These are not the same but that is ok as they treat tau a bit differently.
    # But I think it is more consistent to use calc_iv_masscons:
    warn("Probably better use calc_iv_masscons")

    ud = 2A./(n+1) .* (tau_mean_n./hs_n).^n.*hs_n.^(n+1)
    out = ud.*(1./(1- fsl_n))
    return out*year2sec # m/a
end

"""
$(SIGNATURES)

Calculates ice surface velocity using mass conservation.

Returns iv in [m/a]
"""
function calc_iv_masscons(gl, gb, fsl, qtot_n, hs_n, pp, pm)
    @unpack n, r = pp
    @unpack ele, ws = gb
    fsl_n = on_nodes(fsl)
    ws_n = on_nodes(ws)

    out = similar(hs_n)
    for i=1:length(gb.bands)+1
        if qtot_n[i]<=0
            out[i] = 0
        else
            ubar = qtot_n[i] / (ws_n[i]*hs_n[i])
            out[i] = ubar / ((1-r) * fsl_n[i] + r) *year2sec # m/a
        end
    end
    return out
end

"""
$(SIGNATURES)

Calculates the thickness of one glacier in all elevation bands.
"""
function calc_thick(gb::Bands, qtot_n, fsl, temp, pp::Phys, pm::MPara, pn::Num)
    @unpack relaxfac, relaxmin, niter, reltol, minthick, verbose = pn
    @unpack bands, ws, ele, malphas = gb
    @unpack h_at_top, use_flux_corrector = pm
    @unpack isicecap = gb.gl

    nb = length(bands) # number of cells
    nn = nb + 1 # number of nodes

    fsl_n = on_nodes(fsl)
    temp_n = on_nodes(temp)
    ws_n = on_nodes(ws)

    # Deformation flux normalized by glacier width
    qd_n = flux_deform(qtot_n, fsl_n, pp.r)./ws_n

    # This is some Huss magic:
    if use_flux_corrector
        for i=1:nn
            qd_n[i] = flux_corrector(qd_n[i],i,nn)
        end
    end

    hs0_n = zeros(F,nn)+500  # hs of last step. Fill with IC
    if isnan(h_at_top)
        h_at_top = isicecap ? -1.0 : 0.0
    end
    if h_at_top>=0
        hs0_n[1] = h_at_top
    end
    hs_n = copy(hs0_n)   # hs of current step
    tmp_n = zeros(F,nn)      # work array
    # local tau/h
    # Note: no need to have this on the nodes as an average will be used anyway.
    tauh_local = F[tauh(on_cells(hs_n), gb, ib, pp, pm) for ib=1:nb]
    @assert all(tauh_local.>=0) "tauh_local<0"

    if pn.plotyes
        eval(:(import PyPlot))
        scale = 1e3
        PyPlot.figure()
        PyPlot.subplot(2,1,1)
        PyPlot.plot(gb.xmid/scale,gb.ele, lw=3)
        PyPlot.xlabel("χ (km)")
        PyPlot.ylabel("z (m)")
        PyPlot.hold(true)
        PyPlot.plot(gb.xmid/scale,gb.ele-hs) # IC
        PyPlot.pause(0.05)
        PyPlot.subplot(2,1,2)
        PyPlot.xlabel("χ (km)")
        PyPlot.ylabel("tau (kPa)")
        PyPlot.plot(gb.xmid/scale, tauh_local.*hs/1e3)
        PyPlot.hold(true)
    end
    j=0
    @inbounds @fastmath for j=1:niter
        relax = min(1.0, (1-relaxfac)^j + relaxmin) # relaxation increases with j
        for i_n=2:nn
            hs_n[i_n] = max((thick(hs0_n, qd_n, temp_n, tauh_local, gb, i_n, pp, pm)-hs0_n[i_n])*relax + hs0_n[i_n],
                            minthick)
        end
        if h_at_top>=0 # fixed thickness
            hs0_n[1] = h_at_top
        elseif h_at_top==-1 # keep thickness constant
            hs0_n[1] = hs0_n[2]
        else
            error()
        end
        if any(isnan(hs_n)); error("Got NaN at indices: $(find(isnan(hs_n))) on iteration $j") end
        if norm((hs_n[2:end]-hs0_n[2:end])./hs0_n[2:end],Inf)<reltol; break end
        # update local tau
        for ib=1:nb
            tauh_local[ib] = tauh(on_cells(hs_n,ib), malphas[ib], ws[ib], pp, pm)
            @assert tauh_local[ib]>0 "tauh_local<=0!"
        end
        if pn.plotyes
            PyPlot.subplot(2,1,1)
            PyPlot.plot(gb.xmid/scale,gb.ele-hs)
            PyPlot.subplot(2,1,2)
            PyPlot.plot(gb.xmid/scale, tauh_local.*hs/1e3)
            PyPlot.pause(0.05)
        end

        # swaps bindings
        #        hs, hs0 = hs0, hs
        copy!(tmp_n, hs_n)
        copy!(hs_n, hs0_n)
        copy!(hs0_n, tmp_n)
    end
    if verbose
        if j==niter
            warn("Not converged after $(j) iterations.  Max rel. error of thickness is $(norm((hs_n[2:end]-hs0_n[2:end])./hs0_n[2:end],Inf))")
        else
            println("Converged after $(j) iterations")
        end
    end
    if pn.plotyes
        PyPlot.subplot(2,1,1)
        PyPlot.plot(gb.xmid/scale,gb.ele-hs)
        PyPlot.subplot(2,1,2)
        PyPlot.plot(gb.xmid/scale, tauh_local.*hs/1e3)
    end
    return hs_n, tauh_local.*on_cells(hs_n), j
end


"""
$(SIGNATURES)

Ice thickness at one location for fixed values.

Inputs (on nodes):
- hs_n -- guess for thickness
- qd_n -- deformational ice flux normalized by width
- temp_n -- ice temperature
- tauh_local -- local τ/h (using the previous iteration's h)
- gb -- elevation bands
- i_n -- index of elevation band to calculate
- pp, pm -- parameters

TODO:
- rate factor
- think about staggered vs non-staggered grid (see Clarke&al 2013)
"""
function thick(hs_n::AbstractVector,
               qd_n::AbstractVector,
               temp_n::AbstractVector,
               tauh_local::AbstractVector,
               gb::Bands, i_n::Int,
               pp::Phys, pm::MPara)
    @unpack n = pp
    @unpack mean_tau_dist = pm

    # Eq. 2 of H&F2014:
    if qd_n[i_n]>0
        # tau averaged over suitably many elevation bands:
        tauh_n = bands_mean_fn(tauh_local, hs_n, gb, i_n, mean_tau_dist)
        @assert tauh_n>=0 "tauh<0! at node $i_n"

        A = ice_flow_factor(temp_n[i_n])
        return (qd_n[i_n]*(n+2)
                /
                (2*A*tauh_n^n)
                )^(1/(n+2))
    else  # return zero ice thickness in places where there is
          # negative flux
        return zero(tauh_local[1])
    end
end


# Ice flux
##########


"""
$(SIGNATURES)

Apparent mass balance in meters (elevation) per year [m/a]
"""
apparent_mass_balance(bdot, dhdt) = bdot - dhdt

"""
    on_nodes(band_quantity, [i_n])

Moves a quantity defined on cells (elevation bands) to the nodes
(border points between two bands) by setting:
- node_quantity[1]=band_quantity[1], node_quantity[end]=band_quantity[end]
- averaging for the others
"""
function on_nodes(band_quantity)
    out = zeros(F, length(band_quantity)+1)
    out[1] = band_quantity[1]
    out[end] = band_quantity[end]
    for i_n=2:length(band_quantity)
        out[i_n] = (band_quantity[i_n]+band_quantity[i_n-1])/2
    end
    return out
end
function on_nodes(band_quantity, i_n)
    if i_n==1
        return band_quantity[1]
    elseif i_n==length(band_quantity)+1
        return band_quantity[end]
    else
        return (band_quantity[i_n]+band_quantity[i_n-1])/2
    end
end



"""
    on_cells(node_quantity)

Inverse of `on_node`, uses also averaging.

TODO: this is probably not the best scheme to move h on nodes back to
cells.  Better would maybe be to have two elevation bandings.
"""
function on_cells(node_quantity)
    out = zeros(F, length(node_quantity)-1)
    for i=1:length(node_quantity)-1
        out[i] = (node_quantity[i+1]+node_quantity[i])/2
    end
    return out
end
on_cells(node_quantity, i_c) = (node_quantity[i_c+1]+node_quantity[i_c])/2


"""
$(SIGNATURES)

Calculates the flux between elevations bands from the apparent mass
balance.

Input: btilde [m/a]

Output: flux [m^3 of ice/s]

NOTE: coverts the [/a] of btilde into [/s]!
"""
function flux_n(gb::Bands, btilde::AbstractVector)
    @unpack ls, ws = gb
    qtot_n = zeros(length(gb.bands)+1)
    qtot_n[1] = 0
    for i_n in 2:length(qtot_n)
        qtot_n[i_n] = ls[i_n-1]*ws[i_n-1]*btilde[i_n-1]/year2sec + qtot_n[i_n-1] # [m^3/s]
    end
    return qtot_n
end

"""
$(SIGNATURES)

As flux_n but disallows negative flux

NOTE: coverts the [/a] of btilde into [/s]!
"""
function flux_n_nonnegative(gb::Bands, btilde::AbstractVector)
    @unpack ls, ws = gb
    qtot_n = zeros(length(gb.bands)+1)
    qtot_n[1] = 0
    for i_n in 2:length(qtot_n)
        qtot_n[i_n] = ls[i_n-1]*ws[i_n-1]*btilde[i_n-1]/year2sec + qtot_n[i_n-1] # [m^3/s]
        qtot_n[i_n] = qtot_n[i_n]<0 ? 0*qtot_n[i_n] : qtot_n[i_n]
    end
    return qtot_n
end


"""
$(SIGNATURES)

Where r is pp.r.

Ice flux due to deformation: qd (m^3/a)
"""
flux_deform(qtot_n, fsl_n, r) = qtot_n .*(1- fsl_n./( (1-r).*fsl_n + r))
# function flux_deform(gb::Bands, btilde::AbstractVector, fsl::AbstractVector, pp::Phys, pm::MPara)
#     @unpack ele = gb
#     qtot = flux(gb, btilde, pp)  # TODO fix staggered
#     # move flux to band centers
#     qtot = (qtot[1:end-1]+qtot[2:end])/2
#     return flux_deform.(qtot, fsl, ele)
# end

"""
$(SIGNATURES)

Flux-correction for glacier terminus (?).  Reduces flux at terminus
"""
function flux_corrector(qd,ib,nb)
    ex = 1 # 1 is as Matthias has it
    if nb<=30
        if ib==nb
            qd -= qd*0.6*ex
        elseif ib==nb-1
            qd -= qd*0.4*ex
        elseif ib==nb-2
            qd -= qd*0.2*ex
        end
    else
        if ib==nb
            qd -= qd*0.5*ex
        elseif ib==nb-1
            qd -= qd*0.4*ex
        elseif ib==nb-2
            qd -= qd*0.3*ex
        elseif ib==nb-3
            qd -= qd*0.2*ex
        elseif ib==nb-4
            qd -= qd*0.1*ex
        end
    end
    return qd
end

# Ice dynamics
##############

# """
# Ratio: basal sliding speed / surface speed.

# In:
# ele - surface elevation of point in question


# TODO:
# - update to use a general sliding ratio function
# """
# function ratio_sliding_surface(ele, gl::Glacier, pm::MPara)
#     f1, f2 = pm.fsl_para1, pm.fsl_para2
#     e1, e2 = gl.elemedian, gl.elemin
#     if ele>=e1
#         return f1
#     elseif ele<=e2
#         return f2
#     else
#         return f2 + (ele-e2) * (f1-f2)/(e1-e2)
#     end
# end
# ratio_sliding_surface(gb::Bands, pm::MPara) = [ratio_sliding_surface(e,gb.gl,pm) for e in gb.ele]

"""
$(SIGNATURES)

Shape factor.  Turn off with third argument==false
"""
Fs(h,w,on) = on ? w./(2h+w) : one(h) # thickness, width

"""
$(SIGNATURES)

Basal_shear_stress/thickness: τ/h

Also know as basal drag, in C&P denoted by τ_b.

If there is no shape factor, then this is constant for a given surface
geometry.

TODO: Arguably the tan should be used instead of the sin.
"""
tauh(h::Number, alpha, w, pp::Phys, pm::MPara) =
    Fs(h,w,pm.shapeF).*pp.rho.*pp.g.*sin.(alpha) # == Fs*tau_d/h
tauh(hs::AbstractVector, gb::Bands, ib::Int, pp::Phys, pm::MPara) =
    tauh(hs[ib], gb.malphas[ib], gb.ws[ib], pp, pm)

"""
$(SIGNATURES)

Driving stress
"""
tau_d(h,alpha,w,pp::Phys) = pp.rho*pp.g*sin.(alpha).*h

"""
$(SIGNATURES)

Takes the mean of a band-variable over all bands lying within distance
`+/- dist*hs_n[i_n]`.
"""
function mean_over_bands(var, hs_n, gb, i_n, dist)
    hh = dist*hs_n[i_n]
    @assert hh>=0
    il = searchsortedfirst(gb.xmid, gb.x[i_n]-hh)
    iu = searchsortedlast(gb.xmid, gb.x[i_n]+hh)
    if length(il:iu)>0
        return mean(var[il:iu])
    else
        if il>length(gb.xmid)
            return var[iu]
        elseif iu<1
            return var[il]
        else # return closer
            if abs(gb.xmid[il]-gb.x[i_n]) > abs(gb.xmid[iu]-gb.x[i_n])
                return var[iu]
            else
                return var[il]
            end
        end
    end
end

"""
$(SIGNATURES)

Takes the mean of a node-variable over all nodes lying within distance
`+/- dist*hs_n[i_n]`.
"""
function mean_over_bands_n(var_n, hs_n, gb, i_n, dist)
    hh = dist*hs_n[i_n]
    @assert hh>=0
    il = searchsortedfirst(gb.x, gb.x[i_n]-hh)
    iu = searchsortedlast(gb.x, gb.x[i_n]+hh)
    @assert length(il:iu)>0 # this should always be true
    return mean(var_n[il:iu])
end

"""
$(SIGNATURES)

Takes the mean of a band-variable over all bands lying within distance
`+/- dist*h_`, where h_ is itself mean-ed over `+/- dist*hs_n[i_n]`
"""
function meanmean_over_bands(var, hs_n, gb, i_n, dist)
    hh = mean_over_bands_n(hs_n, hs_n, gb, i_n, dist)*dist
    il = searchsortedfirst(gb.xmid, gb.x[i_n]-hh)
    iu = searchsortedlast(gb.xmid, gb.x[i_n]+hh)
    if length(il:iu)>0
        return mean(var[il:iu])
    else
        if il>length(gb.xmid)
            return var[iu]
        elseif iu<1
            return var[il]
        else # return closer
            if abs(gb.xmid[il]-gb.x[i_n]) > abs(gb.xmid[iu]-gb.x[i_n])
                return var[iu]
            else
                return var[il]
            end
        end
    end
    error("")
end

# Select which of above to use in the mean-tau calculation:
# [use local h, uses an averaged h]
# The latter is much more numerically stable.
const bands_mean_fn = [mean_over_bands, meanmean_over_bands][2]

#######
# Procs
#######

"""
    length1d(x, hs_n)

Length of the glacier, can be Inf.
- if shorter than domain, use that, corrected with a sqrt-fit.
- if longer extrapolate, cleverly (uses a sqrt fit).
"""
function length1d(x, hs_n)
    # Last point occupied by the glacier
    fl = findlast(h->h>1, hs_n)

    # glacier does not exist
    if fl==0
        return zero(x[1])
    end

    ## To get a continuous length fit a
    ## parabola y = a sqrt(x0 - x) through the last few points.

    # take the mean over a few upstream points to make an average upstream point.
    meaner = 3
    xe, xe1 = x[fl], mean(x[max(1,fl-meaner):fl-1])
    he, he1 = hs_n[fl], mean(hs_n[max(1,fl-meaner):fl-1])
    slope = (he-he1) / (xe-xe1)
    if slope>0
        # upsloping end
        if fl==length(x)
            return Inf
        else
            # Upsloping end of glacier not reaching the end:
            # in this case just return to location of the last filled node
            return x[fl]
        end
    else
        # do a parabola fit
        x0 = (xe1 * (he/he1)^2 - xe) / ( (he/he1)^2 - 1)
        # If the extrapolation reaches further than the hs_n==1 then use that
        if fl<length(x) && x[fl+1]<x0
            return x[fl+1]
        else
            return x0
        end
    end
end

##########
# Checks
##########

"""
$(SIGNATURES)

Check that flux is correct
"""
# TODO!
function check_flux1d(gl, gb, tau_mean, qtot, hs, pp)
    @unpack n, r = pp
    tol = 0.01
    for i=1:length(gb.bands)
        u_def = 2A/(n+1) * tau_mean[i]^n*hs[i]^(n+1)
        u_slide = TODO
        u_mean = u_slide + u_def*r
        flux = u_mean*hs[i]
        if !(1-tolabs(flux/qtot[i])<1+tol)
            error("flux not right!")
        end
    end
end
