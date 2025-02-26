import numpy as np
import matplotlib.pyplot as plt


def acoth(x):
    """
       Calculate the inverse hyperbolic cotangent: acoth(x) = 0.5 * log((x+1)/(x-1))
       Note: Only defined for |x| > 1.
       """
    return 0.5 * np.log((x + 1) / (x - 1))


def b_positive_discrete(magnitudes, delta=0.05):
    """
        Estimate b-value via the b-positive method for discretized magnitude data.

        Parameters:
            magnitudes (array-like): A sequence of earthquake magnitudes arranged in chronological order.
            delta (float): The discretization parameter (magnitudes are discretized at intervals of 2*delta, with a default delta=0.05 corresponding to a step size of 0.1).

        Returns:
            float: Estimated b-value. If calculation fails, return None.

        Logic description:
            1. First, discretize the magnitude data: round to the nearest multiple of 2*delta.
               For example, when delta=0.05, 2*delta=0.1, so we retain 1 decimal place.
            2. Calculate the differences between adjacent magnitudes.
            3. Let Mc' = 2*delta. Filter the data to select differences greater than or equal to Mc' (considering numerical tolerance).
            4. Calculate the mean of the differences, then use the formula:

                b = [1/(delta * ln(10))] * acoth((mean(diff) - Mc' + delta)/delta)

               where acoth is the inverse hyperbolic cotangent function.
        """
    decimals = int(round(-np.log10(2 * delta)))
    discretized = np.round(magnitudes, decimals=decimals)

    #  Calculate the difference between neighboring magnitudes
    diffs = np.diff(discretized)

    # Set Mc' = 2*delta and use a small tolerance
    Mc_prime = 2 * delta
    tol = 1e-6
    # Filter out differences greater than or equal to Mc'
    pos_diffs = diffs[diffs >= Mc_prime - tol]

    if len(pos_diffs) == 0:
        print("No positive magnitude difference (>= Mc').")
        b_value = 1000
        return None

    mean_diff = np.mean(pos_diffs)
    # Calculate the input parameters of acoth: (mean_diff - Mc' + delta)/delta
    arg = (mean_diff - Mc_prime + delta) / delta
    if arg <= 1:
        print("The parameter of acoth is less than or equal to 1 and cannot be calculated acoth")
        b_value = 1000
        return None

    b_value = 1.0 / (delta * np.log(10)) * acoth(arg)
    return b_value


def b_positive_discrete_bootstrap(magnitudes, delta=0.05, n_boot=1000):
    """
       Estimate the error of the b_positive_discrete estimate using the bootstrap method.

       Parameters:
           magnitudes (array-like): Original magnitude sequence (arranged in chronological order).
           delta (float): Discretization parameter.
           n_boot (int): Number of bootstrap resampling iterations.

       Returns:
           tuple: (b_mean, b_std)
               b_mean - average b-value estimated from bootstrap
               b_std  - standard deviation of b-values, serving as the error estimate
       """
    boot_b_values = []
    n = len(magnitudes)

    for i in range(n_boot):
        # The original directory is resampled and kept in chronological order
        indices = np.random.choice(np.arange(n), size=n, replace=True)
        indices.sort()  # Ensure that the data is still sorted by time after resampling
        sample = magnitudes[indices]
        b_val_sample = b_positive_discrete(sample, delta=delta)
        if b_val_sample is not None:
            boot_b_values.append(b_val_sample)

    boot_b_values = np.array(boot_b_values)
    if len(boot_b_values) == 0:
        print("Bootstrap A valid B-value estimate could not be obtained")
        return None,None

    return np.mean(boot_b_values), np.std(boot_b_values)


def get_bootstrap_b_values(magnitudes, delta=0.05, n_boot=1000):
    boot_b_values = []
    n = len(magnitudes)
    for i in range(n_boot):
        indices = np.random.choice(np.arange(n), size=n, replace=True)
        indices.sort()
        sample = magnitudes[indices]
        b_val_sample = b_positive_discrete(sample, delta=delta)
        if b_val_sample is not None:
            boot_b_values.append(b_val_sample)
    return np.array(boot_b_values)
