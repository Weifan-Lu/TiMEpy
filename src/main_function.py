import numpy as np
from scipy.optimize import minimize
import visualization as vis
import matplotlib.pyplot as plt


def load_txt_data(filename):
    """
    Read txt file, skip the first line header,
    the file is tab-delimited, and return a two-dimensional NumPy array
    """
    return np.loadtxt(filename, delimiter='\t', skiprows=1)

def calculate_tidal_phase(ttide_cut, StressS_cut, t_ref):
    """
    Calculate tidal phase and related information.

    Inputs:
      ttide_cut   : 1D numpy array, time series
      StressS_cut : 1D numpy array, tidal stress sequence
      t_ref       : Earthquake event occurrence time
    Outputs:
      ephase           : Tidal phase (in radians)
      maxtime          : Time corresponding to the local maximum
      mintime1         : Time corresponding to the local minimum on the left
      mintime2         : Time corresponding to the local minimum on the right
      stress_at_t      : Tidal stress value at the earthquake time
      stress_rate_at_t : Tidal stress rate at the earthquake time
    """
    # Step 1: Compute the stress rate
    stress_rate = np.diff(StressS_cut) / np.diff(ttide_cut)
    # Step 2: Find the index of the earthquake event time (closest index)
    i = np.argmin(np.abs(ttide_cut - t_ref))

    # Step 3: Get the tidal stress value and stress rate at the earthquake time
    stress_at_t = StressS_cut[i]
    if i < len(stress_rate):
        stress_rate_at_t = stress_rate[i]
    else:
        stress_rate_at_t = np.nan

    # Initialize output variables
    ephase = None
    maxtime = None
    mintime1 = None
    mintime2 = None

    # Step 4: Determine the trend of the stress rate
    if stress_rate[i] >= 0:
        # Positive stress rate case
        # Find the local maximum
        j = i
        while j < len(stress_rate) and stress_rate[j] >= 0:
            j += 1
        maxtime = ttide_cut[j - 1]  # Time corresponding to the local maximum

        # Find the local minimum on the right
        k = j
        while k < len(stress_rate) and stress_rate[k] < 0:
            k += 1
        mintime2 = ttide_cut[k - 1]  # Time corresponding to the local minimum on the right

        # Find the local minimum on the left
        l = i
        while l > 0 and stress_rate[l] >= 0:
            l -= 1
        mintime1 = ttide_cut[l + 1]  # Time corresponding to the local minimum on the left

        # Compute the tidal phase
        for m in np.arange(0, 180.5, 0.5):
            ephasetime = (maxtime - mintime1) / 180.0 * m + mintime1
            if ttide_cut[i] <= ephasetime:
                ephase = -180 + m
                break
    else:
        # Negative stress rate case
        # Find the local minimum on the right
        j = i
        while j < len(stress_rate) and stress_rate[j] < 0:
            j += 1
        mintime2 = ttide_cut[j - 1]  # Time corresponding to the local minimum on the right

        # Find the local maximum on the left
        k = i
        while k > 0 and stress_rate[k] < 0:
            k -= 1
        maxtime = ttide_cut[k + 1]  # Time corresponding to the local maximum

        # Find the local minimum on the left
        l = k
        while l > 0 and stress_rate[l] >= 0:
            l -= 1
        mintime1 = ttide_cut[l + 1]  # Time corresponding to the local minimum on the left

        # Compute the tidal phase
        for m in np.arange(0, 180.5, 0.5):
            ephasetime = (mintime2 - maxtime) / 180.0 * m + maxtime
            if ttide_cut[i] <= ephasetime:
                ephase = m
                break

    # Convert tidal phase from degrees to radians
    if ephase is not None:
        ephase = np.deg2rad(ephase)
    else:
        ephase = np.nan

    return ephase, maxtime, mintime1, mintime2, stress_at_t, stress_rate_at_t


def fit_periodic_function(phi_degs, prob):
    """
    Fitting R(φ) = a*cos(φ - φ0) using the least squares method
    Where φ, φ0 are both degrees, and finally get:
    - amplitude = a
    - phase_shift = φ0 (degrees)
    - model_func(φ) can return the fitted R(φ)
    Parameters:
    phi_degs : ndarray, angle (degrees), such as [0, 10, 20, ...]
    prob : ndarray, observation value corresponding to the angle (such as probability density)
    Return:
    amplitude, phase_shift, model_func
    """
    # 1. Design matrix G = [cos(φ), sin(φ)]
    G = np.column_stack((
        np.cos(np.deg2rad(phi_degs)),
        np.sin(np.deg2rad(phi_degs))
    ))

    # 2. Least squares solution [A, B]
    #   R(φ) ~ A*cos(φ) + B*sin(φ)
    p, residuals, rank, s = np.linalg.lstsq(G, prob, rcond=None)
    A, B = p

    # 3. Convert to complex number c = A - iB, so that it is easy to find amplitude and phase at the same time
    c = A - B*1j
    amplitude = np.abs(c)
    phase_shift = -np.angle(c, deg=True)

    # 4. Construct prediction function
    def model_func(phi):
        # Input phi (degrees), output fitted R(φ)
        phi_rad = np.deg2rad(phi)
        return A*np.cos(phi_rad) + B*np.sin(phi_rad)

    return amplitude, phase_shift, model_func

def modulation_phase(data, bin, PhStress, plot_option=None):
    """
    Corresponding to MATLAB's ModulationPhase function

    Input:
    data: phase data (unit: degree)
    bin: the center value of the histogram bin (consistent with the second parameter of hist in MATLAB, requiring uniform distribution)
    PhStress: reference phase data (unit: degree)
    gy: drawing options, such as 'yes2' or 'yes1' (optional)

    Output:
    PM_ra: normalized modulation amplitude
    ph_shift: phase shift (degrees)
    Po: uniform distribution probability density (per bin)
    prob_1: the value of the original histogram after scaling
    """
    data = np.asarray(data)
    bin = np.asarray(bin)
    PhStress = np.asarray(PhStress)

    # MATLAB: [cout,ph1] = hist(data,bin);
    # To be consistent with MATLAB's hist, we need to construct the bin boundaries based on the given bin (bin center):
    if len(bin) < 2:
        raise ValueError("The bin array needs to contain at least two elements to calculate the bin width.")
    bin_edges = np.empty(len(bin) + 1)
    bin_edges[1:-1] = (bin[:-1] + bin[1:]) / 2
    bin_edges[0] = bin[0] - (bin[1] - bin[0]) / 2
    bin_edges[-1] = bin[-1] + (bin[-1] - bin[-2]) / 2
    counts, cc1 = np.histogram(data, bins=bin_edges)
    ph1 = bin.copy()

    diffs = np.diff(bin)
    unique_diffs = np.unique(diffs)
    if unique_diffs.size == 1:
        wbin = unique_diffs[0]
    else:
        raise ValueError("Uneven spacing of bin centers is not supported.")

    Prob_o = counts / (len(data) * wbin)

    counts2, _ = np.histogram(PhStress, bins=bin_edges)

    Po = counts2 / (len(PhStress) * wbin)

    G = np.column_stack((np.cos(np.deg2rad(ph1)), np.sin(np.deg2rad(ph1))))

    Prob = Prob_o / Po - 1

    phi_degs = (bin_edges[:-1] + bin_edges[1:]) / 2
    a_fit, phi0_fit, model = fit_periodic_function(phi_degs ,Prob)
    print(f"Fitted amplitude = {a_fit:.3f}")
    print(f"Fitted phase_shift = {phi0_fit:.3f} deg")
    phi_dense = np.linspace(-180, 180, 30)
    ph_shift  = phi0_fit
    PM_ra = a_fit

    if plot_option is not None and plot_option.lower() == 'yes1':
        vis.plot_phase_modulation(ph1, Prob_o, Po, wbin, PM_ra, ph_shift,model,phi_dense)

    return PM_ra, ph_shift, Po, Prob




def modulation_stress(Stress_AM_bk, Stress_AM, initial_guess, plot_option):
    """
        Corresponding to MATLAB's ModulationStress function
        Stress_AM_S_bk: Reference stress data (one-dimensional array)
        Stress_AM_S: Observed stress data (one-dimensional array)
        initial_guess: Optimize initial parameters, such as [a0, C0]
        plot_option: 'yes' plot, 'no' no plot
        Return: a_estimated, C_estimated, delta_a, delta_c
    """

    # 1. Binning
    edges = np.linspace(np.min(Stress_AM), np.max(Stress_AM), 20)
    counts, _ = np.histogram(Stress_AM, bins=edges)
    counts_bk, _ = np.histogram(Stress_AM_bk, bins=edges)
    event_rate = counts / counts_bk

    # 2. Filter
    nonzero_idx = event_rate>0
    shear_stress = (edges[:-1] + edges[1:]) / 2
    shear_stress = shear_stress[nonzero_idx]

    # 3. Convert to kPa
    shear_stress_kPa = shear_stress / 1000.0

    counts = counts[nonzero_idx]
    counts_bk = counts_bk[nonzero_idx]
    event_rate = counts / counts_bk

    def neg_log_likelihood(params):
        a, C = params
        # To avoid log(0) or negative values, a lower limit is usually added
        model_rate = C * np.exp(a * shear_stress_kPa)
        model_rate = np.where(model_rate <= 0, 1e-60, model_rate)
        return -np.sum(counts * np.log(model_rate) - counts_bk * model_rate)

    def grad_neg_log_likelihood(params):
        a, C = params
        model_rate = C * np.exp(a * shear_stress_kPa)
        model_rate = np.where(model_rate <= 0, 1e-60, model_rate)
        # d(nll)/da
        grad_a = np.sum(
            shear_stress_kPa * (counts_bk * model_rate - counts)
        )
        # d(nll)/dC
        grad_C = np.sum(
            counts_bk * np.exp(a * shear_stress_kPa) - counts / C
        )
        return np.array([grad_a, grad_C])

    def hess_neg_log_likelihood(params):
        a, C = params

        exp_term = np.exp(a * shear_stress_kPa)  # e^{a x_i}
        model_rate = C * exp_term  # C e^{a x_i}

        # H[0,0] = d^2(nll)/(da^2)
        h_aa = np.sum(shear_stress_kPa ** 2 * counts_bk * model_rate)

        # H[1,1] = d^2(nll)/(dC^2)
        h_CC = np.sum(counts / (C ** 2))

        # H[0,1] = H[1,0] = d^2(nll)/(da dC)
        h_aC = np.sum(shear_stress_kPa * counts_bk * exp_term)

        # Make a 2x2 matrix
        return np.array([
            [h_aa, h_aC],
            [h_aC, h_CC]
        ])


    options = {
        'disp': False,  # Whether to print optimization information
        'maxiter': 400, # Maximum number of iterations
        'xtol': 1e-6,  # Or 'tol', 'gtol', etc. can be used for precision adjustment
    }

    res = minimize(
        fun=neg_log_likelihood,
        x0=initial_guess,
        method='Newton-CG',
        jac=grad_neg_log_likelihood,
        hess=hess_neg_log_likelihood,
        options=options
    )
    a_estimated, C_estimated = res.x

    H = hess_neg_log_likelihood(res.x)
    cov_matrix = np.linalg.inv(H)


    std_error_a = np.sqrt(cov_matrix[0, 0])
    std_error_C = np.sqrt(cov_matrix[1, 1])

    delta_a = 2 * std_error_a
    delta_c = 2 * std_error_C


    print("Estimated a:", a_estimated, "with delta:", delta_a)
    print("Estimated C:", C_estimated, "with delta:", delta_c)
    if plot_option == 'yes':
        fig, ax1, ax2 = vis.plot_modulation_stress(
            Stress_AM_bk, Stress_AM, shear_stress, shear_stress_kPa,
            event_rate, a_estimated, C_estimated, delta_a, delta_c
        )


    return a_estimated, C_estimated, delta_a, delta_c, Stress_AM_bk, Stress_AM, shear_stress, shear_stress_kPa, event_rate



def schuster_test(pha):
    """P-value calculation function for Schuster test
        Parameters:
        pha : ndarray
        Phase array in radians

        Returns:
        p_value : float
        P-value of Schuster test
    """
    pha_deg = np.degrees(pha)

    R_x = np.sum(np.cos(np.radians(pha_deg)))
    R_y = np.sum(np.sin(np.radians(pha_deg)))

    R = np.sqrt(R_x ** 2 + R_y ** 2)
    N = len(pha)
    Z = (R ** 2) / N

    p_value = np.exp(-Z)
    return p_value



def grid_modulation_results_2D(Allpha_obs, Allpha_ref, bin_phase, initial_guess, grid_spacing, threshold):
    """
    Performs 2D gridding of the data where each grid cell has a fixed size of grid_spacing degrees.

    Parameters:
        Allpha_obs (ndarray): Observation data with columns representing [?, phase, stress, ..., lat, lon, ...].
                             Assumes that:
                             - Column 2 (index 1) is phase.
                             - Column 3 (index 2) is stress.
                             - Column 11 (index 10) is latitude.
                             - Column 12 (index 11) is longitude.
        Allpha_ref (ndarray): Reference data with similar formatting. Uses:
                             - Column 2 (index 1) for phase.
                             - Column 3 (index 2) for stress.
        bin_phase: Parameter for the ModulationPhase function.
        initial_guess: Initial guess for the ModulationStress function.
        grid_spacing (float): Grid size in degrees (e.g., 0.2).
        threshold (int): Minimum number of samples required in a grid cell to perform calculations.

    Returns:
        grid_centers (ndarray): An array of shape (n, 7) where each row contains:
            [lon_center, lat_center, a_val, p_val, pm_val, pm_pha, number_eq]
    """
    # Extract longitude and latitude from observation data
    # Note: Adjust indices to match your data structure (MATLAB columns 12 and 11 correspond to Python indices 11 and 10)
    lon = Allpha_obs[:, 11]
    lat = Allpha_obs[:, 10]

    # Define grid edges in longitude and latitude directions
    lon_edges = np.arange(np.floor(np.min(lon)), np.ceil(np.max(lon)) + grid_spacing, grid_spacing)
    lat_edges = np.arange(np.floor(np.min(lat)), np.ceil(np.max(lat)) + grid_spacing, grid_spacing)

    # Calculate grid centers
    lon_centers = (lon_edges[:-1] + lon_edges[1:]) / 2
    lat_centers = (lat_edges[:-1] + lat_edges[1:]) / 2

    # Determine the number of grid cells in each direction
    num_bins_lon = len(lon_centers)
    num_bins_lat = len(lat_centers)

    # List to collect grid center results; each row: [lon, lat, a, p, pm, pm_pha, count]
    grid_centers_vol = []

    # Loop over grid cells (lon index first, then lat index)
    for i in range(num_bins_lon):
        for j in range(num_bins_lat):
            # Identify the samples within the current grid cell
            idx = ((lon >= lon_edges[i]) & (lon < lon_edges[i + 1]) &
                   (lat >= lat_edges[j]) & (lat < lat_edges[j + 1]))
            count = np.sum(idx)

            # Skip grid cells with insufficient samples
            if count < threshold:
                continue

            # Extract phase and stress data from the observation data.
            # In MATLAB, column 2 is phase and column 3 is stress.
            phase_obs = Allpha_obs[idx, 1]
            stress_obs = Allpha_obs[idx, 2]

            # Calculate modulation parameters.
            # Convert phases from radians to degrees.
            # Note: The string 'no' is passed as an option (presumably to suppress plotting or messages).

            pm_val, pm_pha, _, _ = modulation_phase(phase_obs * 180 / np.pi,bin_phase,Allpha_ref[:, 1] * 180 / np.pi, 'no')

            # Calculate stress modulation parameter 'a'
            a_val,_, _,_,_, _,_ ,_, _ = modulation_stress(Allpha_ref[:, 2], stress_obs, initial_guess, 'no')

            # Calculate p-value using the Schuster Test
            p_val = schuster_test(phase_obs)

            # Determine the grid cell center coordinates
            lon_center = lon_centers[i]
            lat_center = lat_centers[j]

            # Append the results to the list
            grid_centers_vol.append([lon_center, lat_center, a_val, p_val, pm_val, pm_pha, count])

    # Convert the result to a NumPy array and return
    grid_centers = np.array(grid_centers_vol)
    return grid_centers
