def a_table(Re_S):
    # Define Reynolds number breakpoints and corresponding a1, a2 values
    Re_ranges = [1e1, 1e2, 1e3, 1e4, 1e5]
    a1_vals = [1.400, 1.360, 0.593, 0.321, 0.321]
    a2_vals = [-0.667, -0.657, -0.477, -0.388, -0.388]

    a3 = 1.450  # constant for all Re
    a4 = 0.519  # constant for all Re

    # Determine range bin
    if Re_S < 10:
        idx = 0
    elif Re_S < 1e2:
        idx = 1
    elif Re_S < 1e3:
        idx = 2
    elif Re_S < 1e4:
        idx = 3
    else:
        idx = 4
    return a1_vals[idx], a2_vals[idx], a3, a4