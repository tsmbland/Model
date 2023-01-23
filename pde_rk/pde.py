import numpy as np
from typing import Optional, Tuple


def pde_rk(dxdt: callable, X0: list, Tmax: float, deltat: float, t_eval: np.ndarray,
           killfunc: Optional[callable] = None,
           stabilitycheck: bool = False, maxstep: Optional[float] = None, rk: bool = True) -> Tuple[
    list, float, list, np.ndarray]:
    """
    Function for solving system of PDEs using adaptive Runge-Kutta method
    Adapted from Hubatsch et al., 2019 (see https://github.com/lhcgeneva/PARmodelling)

    Args:
        dxdt: a function that takes list of 1D arrays (one for each species) corresponding to concentrations over space, and returns a list of gradient arrays
        X0: a list specifying the initial state of the system. Will be used as the input to dxdt on the first time step
        Tmax: timepoint at which to terminate simulation
        deltat: initial timestep (this will be adapted throughout the simulation)
        t_eval: a list of timepoints for which to save the state of the system
        killfunc: an optional kill function that takes the same input as dxdt. Integration will terminate when this function returns True
        stabilitycheck: if True, integration will terminate when the system stabilises (changes by less that 1% per 60 time units). Default is False
        maxstep: maximum time step to tolerate
        rk: if True, the function will use an adaptive Runge-Kutta method. If False, the function will use a basic Euler method with a constant time step

    Returns:
        soln: final solution
        time: final time. Will be close to Tmax but not exact due to finite and adaptable time step
        solns: solutions at times specified by t_eval
        times: times corresponding to saved states in solns. Will be close to times specified in t_eval but not exact due to finite and adaptable time step

    """

    # Adaptive step size parameters
    atol = 0.000001
    rtol = 0.000001

    # 5TH ORDER RK COEFFICIENTS for Dormand-Prince
    a21, a31, a32, a41, a42, a43 = 1 / 5, 3 / 40, 9 / 40, 44 / 45, -56 / 15, 32 / 9
    a51, a52, a53, a54 = 19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729
    a61, a62, a63 = 9017 / 3168, -355 / 33, 46732 / 5247
    a64, a65 = 49 / 176, -5103 / 18656
    a71, a72, a73, a74 = 35 / 384, 0, 500 / 1113, 125 / 192
    a75, a76 = -2187 / 6784, 11 / 84

    b1, b2, b3, b4, b5 = 35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784
    b6, b7 = 11 / 84, 0

    bs1, bs2, bs3, bs4 = 5179 / 57600, 0, 7571 / 16695, 393 / 640
    bs5, bs6, bs7 = -92097 / 339200, 187 / 2100, 1 / 40

    # Set up
    X = X0
    X7 = None
    testsoln = dxdt(X)
    nvars = len(testsoln)
    t_stored = []
    X_stored = [[] for _ in range(nvars)]
    n_stored_times = 0
    if maxstep is None:
        maxstep = 10 * deltat

    # Run
    time = 0
    updated = None
    terminate = False
    while not terminate:
        if time > Tmax:
            terminate = True

        if rk is True:
            # Calculate increments for RK45
            if (time == 0) or not updated:
                X1 = dxdt(X)
            else:
                X1 = X7

            X2 = dxdt([X[i] + deltat * (a21 * X1[i]) for i in range(nvars)])
            X3 = dxdt([X[i] + deltat * (a31 * X1[i] + a32 * X2[i]) for i in range(len(X))])
            X4 = dxdt([X[i] + deltat * (a41 * X1[i] + a42 * X2[i] + a43 * X3[i]) for i in range(nvars)])
            X5 = dxdt([X[i] + deltat * (a51 * X1[i] + a52 * X2[i] + a53 * X3[i] + a54 * X4[i]) for i in range(nvars)])
            X6 = dxdt([X[i] + deltat * (a61 * X1[i] + a62 * X2[i] + a63 * X3[i] + a64 * X4[i] + a65 * X5[i]) for i in
                       range(nvars)])
            X7 = dxdt([X[i] + deltat * (a71 * X1[i] + a73 * X3[i] + a74 * X4[i] + a75 * X5[i] + a76 * X6[i]) for i in
                       range(nvars)])

            # Update concentrations using A1-A6 and P1-P6, coefficient for A7 and P7 is 0.
            Xn_new = [X[i] + deltat * (b1 * X1[i] + b3 * X3[i] + b4 * X4[i] + b5 * X5[i] + b6 * X6[i]) for i in
                      range(nvars)]  # b2/7=0

            # Compute difference between fourth and fifth order
            deltaXnerr = [np.max(np.abs(
                (b1 - bs1) * X1[i] + (b3 - bs3) * X3[i] + (b4 - bs4) * X4[i] + (b5 - bs5) * X5[i] + (b6 - bs6) * X6[
                    i] - bs7 * X7[i])) for i in range(nvars)]  # b7 is zero

            # Get maximum concentrations for An and Pn
            yXn = [np.maximum(np.max(np.abs(Xn_new[i])), np.max(np.abs(X[i]))) for i in range(nvars)]

            # Get error scale, combining relative and absolute error
            scaleXn = [atol + yXn[i] * rtol for i in range(nvars)]

            # Compute total error as norm of maximum errors for each species scaled by the error scale
            errs = [(deltaXnerr[i] / scaleXn[i]) ** 2 for i in range(nvars)]
            totalerror = np.sqrt(np.sum(errs) / nvars)

            # Compute new timestep
            # sometimes see "RuntimeWarning: divide by zero encountered in double_scalars". Need to look into
            dtnew = 0.8 * deltat * np.abs(1 / totalerror) ** (1 / 5)

            # Upper and lower bound for timestep to avoid changing too fast
            if dtnew > maxstep:
                dtnew = maxstep
            elif dtnew < deltat / 5:
                dtnew = deltat / 5

            # Compute max percentage change
            change = np.max([np.max(np.abs(X[i] - Xn_new[i]) / Xn_new[i]) * (60 / dtnew) for i in range(nvars)])

            # Set timestep for next round
            deltat = dtnew

        else:
            totalerror = 0
            step = dxdt(X)
            Xn_new = [X[i] + deltat * step[i] for i in range(nvars)]

        # Accept step if error is on the order of error scale or below
        if totalerror < 1:
            time += deltat
            X = Xn_new
            updated = True

            # Store
            if np.sum(time > t_eval) > n_stored_times:
                t_stored.append(time)
                for i in range(nvars):
                    X_stored[i].append(X[i])
                n_stored_times += 1

            # Kill function
            if killfunc is not None:
                brk = killfunc(X)
                if brk:
                    t_stored.append(time)
                    for i in range(nvars):
                        X_stored[i].append(X[i])
                    break

            # Check stability:
            if stabilitycheck:
                if change < 0.001:
                    t_stored.append(time)
                    for i in range(nvars):
                        X_stored[i].append(X[i])
                    break

        else:
            updated = False

    t_stored_ = np.array(t_stored)
    X_stored_ = [np.array(X_stored[i]) for i in range(nvars)]
    return X, time, X_stored_, t_stored_
