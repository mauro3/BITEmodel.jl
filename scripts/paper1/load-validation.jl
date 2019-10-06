include("../rgi-setup.jl")
using BSON, DataStructures
function load_validation(region_str)
    #region_str ="03"
    region = parse(Int, region_str)

    # dir_output = "../output/validation/production-no-dhdt/"
    # dir_output = "../output/validation/production-dhdt/"
    #dir_output = "../output/validation-alps/production-dhdt/"

    if region in [4, 11]
        sub = "apr-20"
        dir_output = "../output/vali/validation-$(region_str)_fit-sigma/$sub"
    else
        sub = "apr-21"
        dir_output = "../output/vali/validation-$(region_str)_fit-sigma/$sub"
    end

    # calibration glaciers
    fl_fit = fit_target -> joinpath(dir_output, "reg$(region_str)_FIT_$(fit_target)")
    # validation glaciers with prior transfer
    fl_test1 = prior_choice -> joinpath(dir_output, "reg$(region_str)_PRIOR_$(prior_choice)-testgls")
    # test glaciers with prior transfer
    fl_test1_fitgls = prior_choice -> joinpath(dir_output, "reg$(region_str)_PRIOR_$(prior_choice)-fitgls")
    # validation glaciers with prior transfer and fit to IV
    fl_test2 = prior_choice -> joinpath(dir_output, "reg$(region_str)_PRIOR_$(prior_choice)_FIT_iv")

    fit_targets = [BM.FitTarget.h, BM.FitTarget.h_iv, BM.FitTarget.iv, BM.FitTarget.length]
    prior_choice = [:base, :h, :iv, :h_iv] # base means no prior transfer

    out_fit = OrderedDict()
    for ft in fit_targets
        tmp = BSON.load(fl_fit(ft)*".bson")[:out_fit]
        out_fit[Symbol(ft)] = tmp
    end

    out_test = OrderedDict()
    for p in prior_choice
        tmp = BSON.load(fl_test1(p)*".bson")[:out_test]
        out_test[(:length, p)] = tmp
    end
    for p in [:base, :h_iv]
        tmp = BSON.load(fl_test2(p)*".bson")[:out_test]
        out_test[(:iv, p)] = tmp
    end

    out_test_fitgls = OrderedDict()
    for p in prior_choice
        tmp = BSON.load(fl_test1_fitgls(p)*".bson")[:out_test]
        out_test_fitgls[(:length, p)] = tmp
    end


    return out_fit, out_test, out_test_fitgls
end
