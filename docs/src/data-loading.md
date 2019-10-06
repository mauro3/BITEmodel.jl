# Data Loading

```@meta
CurrentModule = GV
```

As usual, data loading is the (almost) most complex part of the model.
It tries to be fairly general, maybe too much so.  (Probably, this
should become its own package at some point.)

All the data needed by the model are encoded in the subtypes of
[`DataKind`](@ref), thus listing both at once:

```julia
abstract DataKind
# glacier outline
abstract OutlineData <: DataKind
# misc model parameters
abstract ParaData <: DataKind
# DEM 2D
abstract DEMData <: DataKind
# ice surface velocity
abstract IVData <: DataKind
# Surface mass balance
abstract BdotData <: DataKind
# dh/dt data
abstract DhdtData <: DataKind
# basal sliding fraction
abstract FslData <: DataKind
# ice temperature
abstract TempData <: DataKind
# ice thickness (e.g. radar or nothing)
abstract ThicknessData <: DataKind
# ice flux at terminus
abstract TerminusFluxData <: DataKind
```

Most of these data can be 2D, 1D or 0D.  Except the DEM and outline
need to be 2D, which are (currently) regarded as fixed, i.e. as having
no error (wrong but makes life much easier).  All the others
quantities can/should have an error.  Note that whilst some data will
come from measurements or other, external models, some maybe just
generated (for instance ice temperature).

One data-set is then encoded by
[`DataSet{DK<:DataKind, LOADER}`](@ref), where `LOADER` is a symbol
representing the used loader function (or synthesising function).

TODO: think about mask vs outlines and other such overlapping cases.

# Basic operation
The basic operation is as follows.

Take:
- `gid` a glacier [`GlacierID`](@ref), and
- `pl` load parameters [`LoadPara`](@ref)

for example:
```julia
# choosing Unteraar Glacier from the ITMIX dataset:
gid = GV.ITMIXGlacier{:Unteraar}()
pl = GV.LoadPara()
```

Now use [`load_glacier!`](@ref) to load all the data and parameters:
```julia
gl,pp,pm,pn = load_glacier!(dt,pl,pp,pm)
```
which will load a default set of data for a specific glacier.

This is all which is needed to run the forward and inverse model, see
[Forward Model](@ref) and almost all for the [Inverse Model](@ref).
```julia
hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, vol_ratio, gb =
        GV.fwdm(gl, pp, pm, pn)
```

# Customisation of the loading process

There are several places in which the loading process can be
adjusted.  In order of increasing complexity:

## Setting various options with `LoadPara`

During the load process there are several aspects which can be
customized.

[`LoadPara`](@ref) customization:
- `dataset_opts` holds options for the various `load!` methods
  which are called to load the actual data.  This is probably the most
  used customization.
- `data_root` change path to data directory
- `update_cache` if `==true`, then force an update of the cached data

Not settable/used fields of [`LoadPara`](@ref) are:
- `cache_dir` data-derivatives which are expensive to calculate are
  cached in this folder (e.g. glacier-masks)
- `target_proj` not used yet, but envisaged that data is re-projected
  into this projection (TODO)

Example of setting the standard-deviation by hand for a few
measurement types, as well as setting a few setting in the
model-parameters [`MPara`](@ref):
```julia
gid = ITMIXGlacier{:Unteraar}()
pl = GV.LoadPara(update_cache=false)
# customize pl:
pl.data_root = joinpath(pl.data_root,"ITMIX")
for (k,v) in Dict(GV.ThicknessData=>10.0,
                  GV.BdotData=>0.2,
                  GV.DhdtData=>0.2,
                  GV.IVData=>1.0)
    pl.dataset_opts[k][:sigma]=v
end
pl.dataset_opts[GV.ParaData][:pm] = Dict(
    :bandsize => 30,
    :window_dem_smooth=>100,
    :window_width_smooth=>100, # can make a noisy tau
    :window=>100.0)

# Now load/create the glacier and parameters
gl,pp,pm,pn = load_glacier!(gid,pl)
```

## Choosing data-sources manually

Instead of above, different data-sets can be chosen.  Consider:
```julia
dt = make_datatable(gid, pl)
gl,pp,pm,pn = load_glacier!(dt,pl)
```
which is equivalent to `gl,pp,pm,pn = load_glacier!(gid,pl)`.
However, now the [`DataTable`](@ref) `dt` can be manipulated before
being passed to `load_glacier!`, which allows different datasource to
be selected.

The [`make_datatable`](@ref) function creates a standard datatable for
a given glacier, e.g. for an ITMIX glacier it will select all the
available ITMIX data.  If instead we want to use a different $\dot{b}$
(say), then we could do:
```julia
dt = make_datatable(gid, pl)
# swap out the bdot
bdot_synth = GV.DataSet{GV.BdotData, GV.SyntheticBenchLoader}([""],
                                      pl.dataset_opts[GV.BdotData],
                                                    (Date(),Date()))
dt = GV.DataTable(dt, bdot=bdot_synth)
gl,pp,pm,pn = load_glacier!(dt,pl)
```


# Under the hood

To accommodate as many different file types as possible, as well as
accommodating synthetic data, etc. the [`load_glacier!`](@ref) function calls
into many different methods of the [`load!`](@ref) function which has
the following signature:
```julia
load!(ds::DataSet{<:DataKind, :Loader}, gid::GlacierID, more_args...)
```

The dispatch happens with with first one or two arguments: i.e. a
[`DataSet`](@ref) (parameterized on a `DataKind` and a `Loader`),
potentially additionally dispatching on a glacier-ID too.  An example,
returning the parameters structures, looks like so:
```julia
function load!(ds::DataSet{ParaData}, gid::GlacierID,
               pp::Phys=Phys(), pm::MPara=MPara(), pn::Num=Num())
    kwargs = ds.opts
    pp = Phys(pp; get(kwargs, :pp, Dict{Symbol,Any}())...)
    pm = MPara(pm; get(kwargs, :pm, Dict{Symbol,Any}())...)
    pn = Num(pn; get(kwargs, :pn, Dict{Symbol,Any}())...)
    return pp::Phys, pm::MPara, pn::Num
end
```
Of note is how the options are passed to the constructors `Phys`,
etc., i.e. they are the field `opts` of [`DataSet`](@ref).  Which in
terms gets populated from the field `dataset_opts` is in the
`LoadPara` instance which was used to create the dataset.

Steps:
- create [`load!`](@ref) functions
- create [`make_datatable`](@ref) to create a useful datatable
- add a method to [`special_processing!`](@ref), if needed.
