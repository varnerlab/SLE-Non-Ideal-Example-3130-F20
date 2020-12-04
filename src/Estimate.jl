function obj_function_NRTL(pVector, experimental_data_array, experimental_gamma_array)

    # package parameters into a dictionary -
    model_parameters_dictionary = Dict{String,Any}()
    model_parameters_dictionary["⍺"] = pVector[1]
    model_parameters_dictionary["a12"] = pVector[2]
    model_parameters_dictionary["a21"] = pVector[3]
    model_parameters_dictionary["b12"] = pVector[4]
    model_parameters_dictionary["b21"] = pVector[5]

    # compute NRTL model -
    ga_model = compute_NRTL_gamma_array(experimental_data_array, model_parameters_dictionary)

    # compute the error -
    ga_experimental = experimental_gamma_array[:,1];
    error_vec = (ga_model .- ga_experimental).^2

    # add some "soft" constraints -
    penalty = 100000000*max(0,-1*pVector[1])

    return (sum(error_vec) + penalty)
end

