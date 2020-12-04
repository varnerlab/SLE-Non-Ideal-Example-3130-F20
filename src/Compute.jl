# --
# Given a T value, and single component data, compute the psi value for each component in the mixture.
# --
function compute_psi(T::Float64, data_array::Array{Float64,2})::Array{Float64,1}

    # ok, so we need to compute the psi_* value -
    R = 8.314 # J/mol-K

    # what is the size of the data array?
    # components are on the rows
    # DH = col 1
    # Tm = col 2
    (number_of_components, number_of_cols) = size(data_array)

    # initialize -
    psi_vector = zeros(number_of_components)
    
    # main -
    for component_index = 1:number_of_components
        
        # get data for component_index =
        𝝙H = data_array[component_index,1]
        Tm = data_array[component_index,2]

        # compute psi for component_index -
        psi_value = -(𝝙H/R)*(1/T - 1/Tm)
        #psi_value = (𝝙H/(R*Tm))*((T-Tm)/T)
    
        # package -
        psi_vector[component_index] = exp(psi_value)
    end

    # return -
    return psi_vector
end

# --
# Given an array of T values, and single component data, compute the psi value for each component in the mixture.
# Calls helper function for each value of T
# --
function compute_psi(T_vector::Array{Float64,1},data_array::Array{Float64,2})::Array{Float64,2}

    # initialize -
    number_of_temperature_points = length(T_vector)
    psi_array = Array{Float64,2}(undef,number_of_temperature_points,2)

    # main -
    for (index,T_value) in enumerate(T_vector)
        
        # compute psi_1 and psi_2 -
        psi_vec_tmp = compute_psi(T_value, data_array)
        
        # package -
        psi_array[index,1] = psi_vec_tmp[1]
        psi_array[index,2] = psi_vec_tmp[2]
    end

    # return -
    return psi_array
end

function compute_NRTL_gamma_array(experimemtal_data_array::Array{Float64,2}, 
    NRTL_parameters::Dict{String,Any})::Array{Float64,2}

    # initialize -
    (number_of_temperature_steps, number_of_cols) = size(experimemtal_data_array)
    gamma_array = Array{Float64,2}(undef,number_of_temperature_steps,2)
    
    # get parameters from dictionary -
    a12 = NRTL_parameters["a12"]
    a21 = NRTL_parameters["a21"]
    b12 = NRTL_parameters["b12"]
    b21 = NRTL_parameters["b21"]
    ⍺ = NRTL_parameters["⍺"]

    # main loop -
    for row_index = 1:number_of_temperature_steps

        # get composition -
        x1 = experimemtal_data_array[row_index,2]
        x2 = 1 - x1 # we are binary -

        # get temperature value -
        T = experimemtal_data_array[row_index,1]

        # ok, compute 𝜏12, 𝜏21 etc -
        𝜏12 = a12 + (b12/T)
        𝜏21 = a21 + (b21/T)
        G12 = exp(-1*⍺*𝜏12)
        G21 = exp(-1*⍺*𝜏21)

        # use binary form -
        # gamma 1 -
        term_1 = ((x2)^2)*(𝜏21*(G21/(x1+x2*G21))^2+(𝜏12*G12)/((x2+x1*G12)^2))
        gamma_1 = exp(term_1)

        # gamma 2 -
        term_2 = ((x1)^2)*(𝜏12*(G12/(x2+x1*G12))^2+(𝜏21*G21)/((x1+x2*G21)^2))
        gamma_2 = exp(term_2)

        # package -
        gamma_array[row_index,1] = gamma_1
        gamma_array[row_index,2] = gamma_2
    end

    # return -
    return gamma_array
end

# --
# Computes the composition x1 and z1 given psi value array. Uses SVN eqn 15.20 and 15.21 (derived from SLE matching, ideal phases)
# --
function compute_ideal_composition(psi_array::Array{Float64,2})::Array{Float64,2}

    # initialize -
    (number_of_rows,number_of_cols) = size(psi_array)
    composition_array = Array{Float64,2}(undef,number_of_rows, number_of_cols)
    
    # main -
    # rows = temperatures
    # cols = species 
    for row_index = 1:number_of_rows
        
        # get the 𝜓 values -
        psi_1_value = psi_array[row_index,1]
        psi_2_value = psi_array[row_index,2]

        # compute the compostion -
        # use SVN eqn 15.20 and 15.21 (derived from SLE matching, ideal phases)
        x1 = ((psi_1_value)*(1-psi_2_value))/(psi_1_value - psi_2_value)
        z1 = (1-psi_2_value)/(psi_1_value - psi_2_value)

        # package -
        composition_array[row_index,1] = x1
        composition_array[row_index,2] = z1
    end

    # return -
    return composition_array
end
