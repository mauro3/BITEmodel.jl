"""
Make a symbol with extra bits easily:

   S"dsaf asdfl" => Symbol("dsaf asdfl")
"""
macro S_str(string) Base.Meta.quot(Symbol(string)) end

"""
    @get use_forward_model = opts | 123

is sugar for

    use_forward_model = get(opts, :use_forward_model, 123)
"""
macro get(ex)
    @assert ex.head==:(=)
    @assert ex.args[2].head==:call
    @assert ex.args[2].args[1]==:|
    out = ex.args[1]
    var, def = ex.args[2].args[2:3]
    esc(:($out = get($var, $(QuoteNode(out)), $def)))
end

"""
    @get! use_forward_model = opts | false

is sugar for

    use_forward_model = get!(opts, :use_forward_model, false)
"""
macro get!(ex)
    @assert ex.head==:(=)
    @assert ex.args[2].head==:call
    @assert ex.args[2].args[1]==:|
    out = ex.args[1]
    var, def = ex.args[2].args[2:3]
    esc(:($out = get!($var, $(QuoteNode(out)), $def)))
end

"""
    volume_area_mean_h(gl::Glacier)
    volume_area_mean_h(area::Number; isicecap::Bool=false)
    volume_area_mean_h(area::Number,c::Number,gamma::Number)

Volume Area scaling formula.  Returns mean thickness in meters.
"""
function volume_area_mean_h(area::Number; isicecap::Bool=false)
    if isicecap
        c = 0.034
        gamma = 1.36
    else
        c = 0.054
        gamma = 1.25
    end
    volume_area_mean_h(area,c,gamma)
end
volume_area_mean_h(area::Number,c::Number,gamma::Number) = c*(area/1e6)^(gamma-1) * 1e3
