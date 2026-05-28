
import numpy as np



# ABCD matrix - Boyd Table 4.1
def M_free(d):
    return np.array([[1, d],
                     [0, 1]])

def M_sph_mirror(R):
    return np.array([[1, 0],
                     [-2/R, 1]])


def stability_condition(A,D):
    """ stability condition resonators -1 < (A+D)/2 <1 """
    return -1 <  (A+D)/2   <  1




### Concave Concave cavity 
def concave_concave( roc1, roc2, L_values, wavelength, d= 10e-3, n=1.816   ): 
    """
    Parameters
        roc1, roc2 : float - Radius of curvature of mirror 1 and mirror 2
        L_values : np.array - Array of cavity lengths to test
        wavelength : float 
        d : float - Crystal length in meters. Can be set to 0
        n : float Refractive index of the kktp crystal
    Returns
        np.array(waist_values) : np.array - Array of beam waist values for stable cavity lengths
        np.array(z_waist_values) : np.array - Array of z positions of the waist for stable cavity length
        stability_values : np.array(dtype=bool)  Boolean array indicating cavity stability
    """

    stability_values = np.zeros(L_values.shape, dtype=bool)
    waist_values = []
    z_waist_values = []
    
    for i, L in enumerate(L_values):

        # Propagation matrix 
        L_air = L - d    
        L1 = L_air / 2
        L2 = L_air / 2

        M =  ( M_sph_mirror(roc1) @ M_free(L1) @ M_free(d / n) @ M_free(L2) @ M_sph_mirror(roc2) @ M_free(L2) @ M_free(d / n) @ M_free(L1) )    # starting from the mirror
        #M =  M_free(d / n) @ M_free(L1) @ M_sph_mirror(roc1) @ M_free(L1) @ M_free(d / n) @ M_free(L2) @ M_sph_mirror(roc2) @ M_free(L2) @ M_free(d / n)       # starting from the cristal

        A = M[0, 0]
        B = M[0, 1]
        C = M[1, 0]
        D = M[1, 1]

        # Stability and if, waist
        stability_values[i] = stability_condition(A,D)

        if stability_values[i]:
    
            q_roots = np.roots([C, D - A, -B])
            q0 = q_roots[np.imag(q_roots) > 0][0]

            z_array = np.linspace(0, L, 1000)
            q_out = q0 + z_array  
            
            w_array = np.sqrt(wavelength / (np.pi * np.imag( -1 / q_out)))

            i0 = np.argmin(w_array)

        
            waist_values.append(w_array[i0] )   # m
            z_waist_values.append(z_array[i0] )    # m 
            
    return np.array(waist_values), np.array(z_waist_values),  stability_values



# Kogelnik, Laser beams and resonators, Table II
def kogelnik_stability_values( R1, R2, L_values):
    stability_values = np.zeros(L_values.shape, dtype=bool)
    
    for i, L in enumerate(L_values):

        propagation_matrix =  M_sph_mirror(R1) @ M_free(L) @ M_sph_mirror(R2)  @ M_free(L)
        A = propagation_matrix[0, 0]
        B = propagation_matrix[0, 1]
        C = propagation_matrix[1, 0]
        D = propagation_matrix[1, 1]
        stability_values[i] = stability_condition(A,D)

    return stability_values


def kogelnik_waist_concave_concave(R1, R2, L, wavelength): 

    x =  np.sqrt(L*(R1-L)*(R2-L)*(R1+R2-L))  / np.abs(R1+R2 -2*L)  # = b/2
    waist_carre = wavelength * x / np.pi 
    waist = np.sqrt( waist_carre )
    return waist




### Plan Concave
def plan_concave( roc1, L_values, wavelength, wz_target, beam_profile = False): 

    n = 1.816       # Indice optique cristal
    d = 10e-3       # cristal size 10 mm

    stability_values = np.zeros(L_values.shape, dtype=bool)

    waist_values = []
    wz_close_to_taget_values = []
    z_close_to_taget_values = []
    i_target = None
    
    wz__lenths_array = []
    z__lenths_array = []
    


    for i, L in enumerate(L_values):

        # Propagation matrix 
        M = ( M_sph_mirror(roc1) @ M_free(L) @ M_free(L))
        A = M[0, 0]
        B = M[0, 1]
        C = M[1, 0]
        D = M[1, 1]

        # Stability and if, waist
        stability_values[i] = stability_condition(A,D)
        z_array = np.linspace(0, L, 100)

        if stability_values[i]:
            q_roots = np.roots([C, D - A, -B])
            q0 = q_roots[np.imag(q_roots) > 0][0]

            w0 = np.sqrt(wavelength / (np.pi * np.imag(-1 / q0))) # waist is at z=0 (plane mirror) by definition
            waist_values.append(w0 * 1e6)

            # w(z) target value given in the fonction in m
            q_out = q0 + z_array  
            w_array = np.sqrt(wavelength / (np.pi * np.imag( -1 / q_out)))

            i_target = np.argmin(np.abs(w_array - wz_target))
            wz_close_to_taget_values.append( w_array[i_target])
            z_close_to_taget_values.append( z_array[i_target])

            if beam_profile: 
                wz__lenths_array.append(np.array(w_array))
                z__lenths_array.append(np.array(z_array))


    return np.array(waist_values),    stability_values,    np.array(wz_close_to_taget_values),      np.array(z_close_to_taget_values),      wz__lenths_array,    z__lenths_array