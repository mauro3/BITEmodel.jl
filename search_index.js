var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "BITE-Model: Bayesian Ice Thickness Estimation Model",
    "title": "BITE-Model: Bayesian Ice Thickness Estimation Model",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#BITE-Model:-Bayesian-Ice-Thickness-Estimation-Model-1",
    "page": "BITE-Model: Bayesian Ice Thickness Estimation Model",
    "title": "BITE-Model: Bayesian Ice Thickness Estimation Model",
    "category": "section",
    "text": "(Image: Vernagtferner)Surging, biting Vernagtferner by R. Reschreiter (1911)This model can be used to estimate ice thickness maps of mountain glaciers and ice caps.  It needs surface elevation and a glacier outline and can make use of the following additional data: surface mass balance, elevation change, ice thickness measurements (usually radar), and also expert guesses.  It uses a Bayesian framework implemented numerically using a Markov chain Monte Carlo method."
},

{
    "location": "index.html#TOC-1",
    "page": "BITE-Model: Bayesian Ice Thickness Estimation Model",
    "title": "TOC",
    "category": "section",
    "text": "Pages = [\n    \"data-loading.md\",\n    \"forward.md\",\n    \"inverse.md\",\n    \"running.md\",\n    \"api.md\",\n    \"index-index.md\",\n]\nDepth = 1"
},

{
    "location": "data-loading.html#",
    "page": "Data Loading",
    "title": "Data Loading",
    "category": "page",
    "text": ""
},

{
    "location": "data-loading.html#Data-Loading-1",
    "page": "Data Loading",
    "title": "Data Loading",
    "category": "section",
    "text": "CurrentModule = GVAs usual, data loading is the (almost) most complex part of the model. It tries to be fairly general, maybe too much so.  (Probably, this should become its own package at some point.)All the data needed by the model are encoded in the subtypes of DataKind, thus listing both at once:abstract DataKind\n# glacier outline\nabstract OutlineData <: DataKind\n# misc model parameters\nabstract ParaData <: DataKind\n# DEM 2D\nabstract DEMData <: DataKind\n# ice surface velocity\nabstract IVData <: DataKind\n# Surface mass balance\nabstract BdotData <: DataKind\n# dh/dt data\nabstract DhdtData <: DataKind\n# basal sliding fraction\nabstract FslData <: DataKind\n# ice temperature\nabstract TempData <: DataKind\n# ice thickness (e.g. radar or nothing)\nabstract ThicknessData <: DataKind\n# ice flux at terminus\nabstract TerminusFluxData <: DataKindMost of these data can be 2D, 1D or 0D.  Except the DEM and outline need to be 2D, which are (currently) regarded as fixed, i.e. as having no error (wrong but makes life much easier).  All the others quantities can/should have an error.  Note that whilst some data will come from measurements or other, external models, some maybe just generated (for instance ice temperature).One data-set is then encoded by DataSet{DK<:DataKind, LOADER}, where LOADER is a symbol representing the used loader function (or synthesising function).TODO: think about mask vs outlines and other such overlapping cases."
},

{
    "location": "data-loading.html#Basic-operation-1",
    "page": "Data Loading",
    "title": "Basic operation",
    "category": "section",
    "text": "The basic operation is as follows.Take:gid a glacier GlacierID, and\npl load parameters LoadParafor example:# choosing Unteraar Glacier from the ITMIX dataset:\ngid = GV.ITMIXGlacier{:Unteraar}()\npl = GV.LoadPara()Now use load_glacier! to load all the data and parameters:gl,pp,pm,pn = load_glacier!(dt,pl,pp,pm)which will load a default set of data for a specific glacier.This is all which is needed to run the forward and inverse model, see Forward Model and almost all for the Inverse Model.hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, vol_ratio, gb =\n        GV.fwdm(gl, pp, pm, pn)"
},

{
    "location": "data-loading.html#Customisation-of-the-loading-process-1",
    "page": "Data Loading",
    "title": "Customisation of the loading process",
    "category": "section",
    "text": "There are several places in which the loading process can be adjusted.  In order of increasing complexity:"
},

{
    "location": "data-loading.html#Setting-various-options-with-LoadPara-1",
    "page": "Data Loading",
    "title": "Setting various options with LoadPara",
    "category": "section",
    "text": "During the load process there are several aspects which can be customized.LoadPara customization:dataset_opts holds options for the various load! methods which are called to load the actual data.  This is probably the most used customization.\ndata_root change path to data directory\nupdate_cache if ==true, then force an update of the cached dataNot settable/used fields of LoadPara are:cache_dir data-derivatives which are expensive to calculate are cached in this folder (e.g. glacier-masks)\ntarget_proj not used yet, but envisaged that data is re-projected into this projection (TODO)Example of setting the standard-deviation by hand for a few measurement types, as well as setting a few setting in the model-parameters MPara:gid = ITMIXGlacier{:Unteraar}()\npl = GV.LoadPara(update_cache=false)\n# customize pl:\npl.data_root = joinpath(pl.data_root,\"ITMIX\")\nfor (k,v) in Dict(GV.ThicknessData=>10.0,\n                  GV.BdotData=>0.2,\n                  GV.DhdtData=>0.2,\n                  GV.IVData=>1.0)\n    pl.dataset_opts[k][:sigma]=v\nend\npl.dataset_opts[GV.ParaData][:pm] = Dict(\n    :bandsize => 30,\n    :window_dem_smooth=>100,\n    :window_width_smooth=>100, # can make a noisy tau\n    :window=>100.0)\n\n# Now load/create the glacier and parameters\ngl,pp,pm,pn = load_glacier!(gid,pl)"
},

{
    "location": "data-loading.html#Choosing-data-sources-manually-1",
    "page": "Data Loading",
    "title": "Choosing data-sources manually",
    "category": "section",
    "text": "Instead of above, different data-sets can be chosen.  Consider:dt = make_datatable(gid, pl)\ngl,pp,pm,pn = load_glacier!(dt,pl)which is equivalent to gl,pp,pm,pn = load_glacier!(gid,pl). However, now the DataTable dt can be manipulated before being passed to load_glacier!, which allows different datasource to be selected.The make_datatable function creates a standard datatable for a given glacier, e.g. for an ITMIX glacier it will select all the available ITMIX data.  If instead we want to use a different dotb (say), then we could do:dt = make_datatable(gid, pl)\n# swap out the bdot\nbdot_synth = GV.DataSet{GV.BdotData, GV.SyntheticBenchLoader}([\"\"],\n                                      pl.dataset_opts[GV.BdotData],\n                                                    (Date(),Date()))\ndt = GV.DataTable(dt, bdot=bdot_synth)\ngl,pp,pm,pn = load_glacier!(dt,pl)"
},

{
    "location": "data-loading.html#Under-the-hood-1",
    "page": "Data Loading",
    "title": "Under the hood",
    "category": "section",
    "text": "To accommodate as many different file types as possible, as well as accommodating synthetic data, etc. the load_glacier! function calls into many different methods of the load! function which has the following signature:load!(ds::DataSet{<:DataKind, :Loader}, gid::GlacierID, more_args...)The dispatch happens with with first one or two arguments: i.e. a DataSet (parameterized on a DataKind and a Loader), potentially additionally dispatching on a glacier-ID too.  An example, returning the parameters structures, looks like so:function load!(ds::DataSet{ParaData}, gid::GlacierID,\n               pp::Phys=Phys(), pm::MPara=MPara(), pn::Num=Num())\n    kwargs = ds.opts\n    pp = Phys(pp; get(kwargs, :pp, Dict{Symbol,Any}())...)\n    pm = MPara(pm; get(kwargs, :pm, Dict{Symbol,Any}())...)\n    pn = Num(pn; get(kwargs, :pn, Dict{Symbol,Any}())...)\n    return pp::Phys, pm::MPara, pn::Num\nendOf note is how the options are passed to the constructors Phys, etc., i.e. they are the field opts of DataSet.  Which in terms gets populated from the field dataset_opts is in the LoadPara instance which was used to create the dataset.Steps:create load! functions\ncreate make_datatable to create a useful datatable\nadd a method to special_processing!, if needed."
},

{
    "location": "forward.html#",
    "page": "Forward Model",
    "title": "Forward Model",
    "category": "page",
    "text": ""
},

{
    "location": "forward.html#Forward-Model-1",
    "page": "Forward Model",
    "title": "Forward Model",
    "category": "section",
    "text": ""
},

{
    "location": "inverse.html#",
    "page": "Inverse Model",
    "title": "Inverse Model",
    "category": "page",
    "text": ""
},

{
    "location": "inverse.html#Inverse-Model-1",
    "page": "Inverse Model",
    "title": "Inverse Model",
    "category": "section",
    "text": ""
},

{
    "location": "postproc.html#",
    "page": "Post-processing and Plotting",
    "title": "Post-processing and Plotting",
    "category": "page",
    "text": ""
},

{
    "location": "postproc.html#Post-processing-and-Plotting-1",
    "page": "Post-processing and Plotting",
    "title": "Post-processing and Plotting",
    "category": "section",
    "text": ""
},

{
    "location": "api.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": ""
},

{
    "location": "api.html#Data-Loading-1",
    "page": "API",
    "title": "Data Loading",
    "category": "section",
    "text": "Modules = [GV]\nPages = [\"*/data-loading.jl\", \"*/data-loading-misc.jl\"]"
},

{
    "location": "api.html#Loaders-1",
    "page": "API",
    "title": "Loaders",
    "category": "section",
    "text": "Modules = [GV]\nPages = [\"*/loaders/*\"]"
},

{
    "location": "api.html#Forward-model-1",
    "page": "API",
    "title": "Forward model",
    "category": "section",
    "text": "Modules = [GV]\nPages = [\"*/OneD.jl\", \"*/TwoD.jl\"]"
},

{
    "location": "api.html#Inverse-model-1",
    "page": "API",
    "title": "Inverse model",
    "category": "section",
    "text": "Modules = [GV]\nPages = [\"*/mcmc/*\"]"
},

{
    "location": "api.html#Utils-1",
    "page": "API",
    "title": "Utils",
    "category": "section",
    "text": "Modules = [GV]\nPages = [\"*/postproc.jl\", \"*/plotting.jl\", \"*/plotting-inverse.jl\"]"
},

{
    "location": "api.html#Misc-1",
    "page": "API",
    "title": "Misc",
    "category": "section",
    "text": "Modules = [GV]\nPages = [\"*/GV.jl\",\"*/misc.jl\", \"*/data-structs.jl\"]"
},

{
    "location": "index-index.html#",
    "page": "Index",
    "title": "Index",
    "category": "page",
    "text": ""
},

{
    "location": "index-index.html#Index-1",
    "page": "Index",
    "title": "Index",
    "category": "section",
    "text": ""
},

]}
