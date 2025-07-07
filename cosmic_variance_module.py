import numpy as np
import pandas as pd
mean_z = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.5, 1.9, 2.5, 3.5]
delta_z = 0.2
mass_array = [8.75, 9.25, 9.75, 10.25, 10.75, 11.25]
def calc_COSMOS_cosmic_variance(mean_z, delta_z, log_m_stellar, survey="GOODS"):
    """
    Calculates the cosmic variance in the COSMOS field, based on the 'recipe' in:
    "A COSMIC VARIANCE COOKBOOK" by Moster+2011.
    """
    # Step 1: Choose survey field in Table 1 (values below are from Table 3)
    if survey == "UDF":
        sigma_a, sigma_b, beta = 0.251, 0.364, 0.358
    elif survey == "GOODS":
        sigma_a, sigma_b, beta = 0.261, 0.854, 0.684
    elif survey == "GEMS":
        sigma_a, sigma_b, beta = 0.161, 0.520, 0.729
    elif survey == "EGS":
        sigma_a, sigma_b, beta = 0.128, 0.383, 0.673
    elif survey == "COSMOS":
        sigma_a, sigma_b, beta = 0.069, 0.234, 0.834
    else:
        raise ValueError("Survey field must be specified correctly.")

    # Step 4a: Calculate dark matter root cosmic variance
    sigma_DM_meanz_deltaz02 = sigma_a / (mean_z**beta + sigma_b)  # Equation (10)
    sigma_DM = sigma_DM_meanz_deltaz02
    # Step 4b: Calculate the galaxy bias
    dict_galaxy_bias_fit_parameters = {
        8.75: [0.062, 2.59, 1.025],
        9.25: [0.074, 2.58, 1.039],
        9.75: [0.042, 3.17, 1.147],
        10.25: [0.053, 3.07, 1.225],
        10.75: [0.069, 3.19, 1.269],
        11.25: [0.173, 2.89, 1.438],
        '>11.0': [0.185, 2.86, 1.448]
    }

    if log_m_stellar not in dict_galaxy_bias_fit_parameters:
        log_m_stellar = min(dict_galaxy_bias_fit_parameters.keys(), key=lambda x: abs(x - log_m_stellar))

    b_0, b_1, b_2 = dict_galaxy_bias_fit_parameters[log_m_stellar]
    b_Mstellar_meanz = b_0 * (mean_z + 1.0)**b_1 + b_2  # Equation (13)

    # Step 4c: Compute the root cosmic variance for galaxies
    delta_gg = b_Mstellar_meanz * sigma_DM_meanz_deltaz02 * np.sqrt(0.2 / delta_z)

    return delta_gg, sigma_DM

def get_cosmic_variance_array(mean_z, delta_z, mass_array, survey="GOODS"):
    """
    Uses "calc_COSMOS_cosmic_variance" to get the cosmic variance for each combination of
    mean_z and log(stellar mass).
    """
    sigma_DM_List = []
    cosmic_var_matrix = np.zeros((len(mean_z), len(mass_array)))

    for i, z in enumerate(mean_z):
        for j, mass in enumerate(mass_array):
            mod = round(round(mass, 1) % 0.5, 2)
            if mod <= 0.0001:
                mass -= 0.25
                mass = round(mass, 2)
            else:
                mod_ = round(round(mass, 1) % 0.25, 2)
                mass_rounded = (mass // 0.25) * 0.25
                if round(round(mass_rounded, 2) % 0.5, 2) <= 0.0001:
                    mass_rounded += 0.25
                mass = mass_rounded

            if mass > 11.25:
                mass = '>11.0'

            cosmic_var_matrix[i, j],sigma_dm = calc_COSMOS_cosmic_variance(z, delta_z, log_m_stellar=mass, survey=survey)
            sigma_DM_List.append(sigma_dm)
            sigma_DM = sigma_DM_List[:: 6]
    return cosmic_var_matrix, sigma_DM
