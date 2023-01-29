import math as mt
import cmath as cm
import matplotlib.pyplot as plt

pi: float = mt.pi
j: complex = complex(0, 1)
inf: complex = complex("inf")
nan: complex = complex("nan")


def isnan(z: float | complex) -> bool:
    """
    Verifies if the input is NaN (Not a Number). If both the real and imaginary parts of z are not NaN, returns False. Returns True otherwise.

    Args:
        z (float | complex): Number to be tested.

    Returns:
        bool: Boolean result of test.
    """

    a: float = complex(z).real
    b: float = complex(z).imag

    return mt.isnan(a) or mt.isnan(b)


def isinf(z: float | complex) -> bool:
    """
    Verifies if the input is infinite. If z is NaN or both the real and imaginary parts of z are not inifnite, returns False. Returns True otherwise.

    Args:
        z (float | complex): Number to be tested.

    Returns:
        bool: Boolean result of test.
    """

    if isnan(z):
        return False

    a: float = complex(z).real
    b: float = complex(z).imag

    return mt.isinf(a) or mt.isinf(b)


def isfinite(z: float | complex) -> bool:
    """
    Verifies if the input is finite. If z is either infinite or NaN returns False. Returns True otherwise.

    Args:
        z (float | complex): Number to be tested.

    Returns:
        bool: Boolean result of test.
    """

    return not isinf(z) and not isnan(z)


def dB(z: int | float | complex) -> float:
    """
    Calculates the magnitude of z in decibels. If z is null or infinite, returns the float versions of -inf and inf respectively. Returns the usual power ratio decibel formula (10log|z|) otherwise.

    Args:
        z (int | float | complex): Value to be converted into decibels.

    Returns:
        float: Result in decibels.
    """

    if z == 0:
        return -inf.real

    if isinf(z):
        return inf.real

    return 10*mt.log10(abs(z))


def real(z: int | float | complex) -> float:
    """
    Calculate the real part of the input.

    Args:
        z (int | float | complex): Input.

    Returns:
        float: Real part of the input.
    """

    return complex(z).real


def imag(z: int | float | complex) -> float:
    """
    Calculate the imaginary part of the input.

    Args:
        z (int | float | complex): Input.

    Returns:
        float: Imaginary part of the input.
    """

    return complex(z).imag


def conj(z: int | float | complex) -> complex:
    """
    Calculates the complex conjugate of the input.

    Args:
        z (int | float | complex): Input.

    Returns:
        float: Conjugate of the input.
    """

    a: float = real(z)
    b: float = imag(z)

    return complex(a, -b)


def pol(r: int | float, theta: int | float) -> complex:
    """
    Calculates the complex number given it's magnitude and phase in radians.

    Args:
        r (int | float): Magnitude of the number.
        theta (int | float): Phase of the number in radians.

    Returns:
        complex: Resulting complex number.
    """

    return r*cm.exp(j*theta)


def resistance_circle(r: int | float, d: float = 0.1, N_max: int = 100) -> None:
    """
    Plots a circle of constant resistance on the Smith chart. Does not show the plot.

    Args:
        r (int | float): Normalized resistance.
        d (float, optional): Distance between consecutive points. If the distance is such that more than N_max points are needed, limit it at N_max. Defaults to 0.1.
        N_max (int, optional): Maximum amount of points. Defaults to 100.
    """

    if isinf(r):
        plt.plot(1, 0, color="grey")

    R = 1/(1+r)

    N: int = min(N_max, 1+mt.floor(2*pi/mt.acos(1-(d**2)/(2*R**2))))

    theta_list: list[float] = [(2*pi*n)/(N-1) for n in range(N)]

    x_list: list[float] = [r/(r+1) + R*mt.cos(theta) for theta in theta_list]
    y_list: list[float] = [R*mt.sin(theta) for theta in theta_list]

    plt.plot(x_list, y_list, color="grey")


def reactance_circle(x: int | float, d: float = 0.1, N_max: int = 100) -> None:
    """
    Plots a circle of constant reactance on the Smith chart. Does not show the plot.

    Args:
        x (int | float): Normalized reactance.
        d (float, optional): Distance between consecutive points. If the distance is such that more than N_max points are needed, limit it at N_max. Defaults to 0.1.
        N_max (int, optional): Maximum amount of points. Defaults to 100.
    """

    if isinf(x):
        plt.plot(1, 0, color="grey")

    theta_min: float = 0
    theta_max: float = 0

    if x < 0:
        theta_min = 3*pi/2
        theta_max = 2*pi+2*mt.atan((x+1)/(x-1))
    elif x == 0:
        plt.plot((-1, 1), (0, 0), color="grey")
        return
    elif 0 < x < 1:
        theta_min = 2*pi+2*mt.atan((x+1)/(x-1))
        theta_max = 3*pi/2
    elif x == 1:
        theta_min = pi
        theta_max = 3*pi/2
    else:
        theta_min = 2*mt.atan((x+1)/(x-1))
        theta_max = 3*pi/2

    d_theta: float = theta_max - theta_min

    R = 1/x

    N: int = min(N_max, 1+mt.floor(d_theta/mt.acos(1-(d**2)/(2*R**2))))

    theta_list: list[float] = [theta_min+(d_theta*n)/(N-1) for n in range(N)]

    x_list: list[float] = [1 + R*mt.cos(theta) for theta in theta_list]
    y_list: list[float] = [1/x + R*mt.sin(theta) for theta in theta_list]

    plt.plot(x_list, y_list, color="grey")


def plot_smith(d: float = 0.1, N_max: int = 1000) -> None:
    """
    Plots a Smith chart with impedance lines. Does not show the plot.

    Args:
        d (float, optional): Distance between consecutive points. If the distance is such that more than N_max points are needed per circle, limit it at N_max. Defaults to 0.1.
        N_max (int, optional): Maximum amount of points for all circlescombined. Defaults to 1000.
    """

    plt.axis(False)
    plt.axis((-1.7, 1.7, -1.2, 1.2))
    plt.grid(False)

    for r in (0, 1/3, 1, 3):
        resistance_circle(r, d, N_max)

    for x in (-2, -1, 1-mt.sqrt(2), 0, mt.sqrt(2)-1, 1, 2):
        reactance_circle(x, d, N_max)


def z_to_s(Z: int | float | complex, Z0: int | float | complex = 50) -> complex:
    """
    Converts the given impedance to it's corresponding reflection coefficient.

    Args:
        Z (int | float | complex): Impedance to be converted in ohms.
        Z0 (int | float | complex, optional): Reference impedance in ohms. Defaults to 50.

    Returns:
        complex: Corresponding reflection coefficient.
    """

    if Z+Z0 == 0:
        return inf

    if isinf(Z):
        return 1

    return (Z-Z0)/(Z+Z0)


def s_to_z(gamma: int | float | complex, Z0: int | float | complex = 50) -> complex:
    """
    Converts the given reflection coefficient to it's corresponding impedance.

    Args:
        gamma (int | float | complex): Reflection coefficient to be converted.
        Z0 (int | float | complex, optional): Reference impedance in ohms. Defaults to 50.

    Returns:
        complex: Corresponding impedance.
    """

    if gamma == 1:
        return inf

    if isinf(gamma):
        return -Z0

    return Z0*(1+gamma)/(1-gamma)


def input_reflection(S: tuple[complex, complex, complex, complex], gammaL: int | float | complex) -> complex:
    """
    Calculates the reflection coefficient at the input of a 2-port network when the output is loaded.

    Args:
        S (tuple[complex, complex, complex, complex]): Unfolded scattering matrix (S11, S21, S12, S22).
        gammaL (int | float | complex): Reflection coefficient of the impedance loading the output.

    Returns:
        complex: Reflection coefficient at the input.
    """

    # TODO not workin for passtrough networks

    S11: complex = S[0]
    S21: complex = S[1]
    S12: complex = S[2]
    S22: complex = S[3]

    if 1-S22*gammaL == 0 and S21*S12*gammaL != 0:
        return inf

    if 1-S22*gammaL == 0 and S21*S12*gammaL == 0:
        return nan

    return S11 + (S21*S12*gammaL)/(1-S22*gammaL)


def output_reflection(S: tuple[complex, complex, complex, complex], gammaS: int | float | complex) -> complex:
    """
    Calculates the reflection coefficient at the output of a 2-port network when the input is loaded.

    Args:
        S (tuple[complex, complex, complex, complex]): Unfolded scattering matrix (S11, S21, S12, S22).
        gammaS (int | float | complex): Reflection coefficient of the impedance loading the input.

    Returns:
        complex: Reflection coefficient at the output.
    """

    # TODO not workin for passtrough networks

    S11: complex = S[0]
    S21: complex = S[1]
    S12: complex = S[2]
    S22: complex = S[3]

    if 1-S11*gammaS == 0 and 1-S11*gammaS != 0:
        return inf

    if 1-S11*gammaS == 0 and 1-S11*gammaS == 0:
        return nan

    return S22 + (S21*S12*gammaS)/(1-S11*gammaS)


def rollet(S: tuple[complex, complex, complex, complex]) -> float:
    """
    Calculates the Rollet stabilty factor of a 2-port network.

    Args:
        S (tuple[complex, complex, complex, complex]): Unfolded scattering matrix (S11, S21, S12, S22).

    Returns:
        float: Rollet stability factor of the network.
    """

    S11: complex = S[0]
    S21: complex = S[1]
    S12: complex = S[2]
    S22: complex = S[3]
    delta: complex = S11*S22-S12*S21

    numer = 1+abs(delta)**2 - abs(S11)**2 - abs(S22)**2
    denom = 2*abs(S12)*abs(S21)

    if denom == 0 and numer > 0:
        return inf.real

    if denom == 0 and numer < 0:
        return -inf.real

    if denom == 0 and numer == 0:
        return nan.real

    return numer/denom


def mu_stab(S: tuple[complex, complex, complex, complex]) -> tuple[float, float]:
    """
    Calculates the mu stability factors of a 2-port network.

    Args:
        S (tuple[complex, complex, complex, complex]): Unfolded scattering matrix (S11, S21, S12, S22).

    Returns:
        tuple[float, float]: (muS, muL)
            muS (float): Mu factor in source plane.
            muL (float): Mu factor in load plane.
    """

    S11: complex = S[0]
    S21: complex = S[1]
    S12: complex = S[2]
    S22: complex = S[3]
    delta: complex = S11*S22-S12*S21

    numerS: float = 1 - abs(S22)**2
    denomS: float = abs(S11 - conj(S22)*delta) + abs(S21*S12)

    if denomS == 0 and numerS != 0:
        muS = inf.real
    elif denomS == 0 and numerS == 0:
        muS = nan.real
    else:
        muS = numerS/denomS

    numerL: float = 1 - abs(S11)**2
    denomL: float = abs(S22 - conj(S11)*delta) + abs(S21*S12)

    if denomL == 0 and numerL != 0:
        muL = inf.real
    elif denomL == 0 and numerL == 0:
        muL = nan.real
    else:
        muL = numerL/denomL

    return (muS, muL)


def gain_transducer(S: tuple[complex, complex, complex, complex], gammaS: int | float | complex, gammaL: int | float | complex) -> float:
    """
    Calculates the transducer power gain of a 2-port network.

    Args:
        S (tuple[complex, complex, complex, complex]): Unfolded scattering matrix (S11, S21, S12, S22).
        gammaS (int | float | complex): Reflection coefficient of the source.
        gammaL (int | float | complex): Reflection coefficient of the load.

    Returns:
        float: Transducer power gain of the network in linear scale.
    """

    S11: complex = S[0]
    S21: complex = S[1]
    S12: complex = S[2]
    S22: complex = S[3]

    if S21 == 0:
        return 0

    numer: float = (1-abs(gammaS)**2)*(abs(S21)**2)*(1-abs(gammaL)**2)
    denom: float = abs((1-S11*gammaS)*(1-S22*gammaL)-S12*S21*gammaL*gammaS)**2

    if denom == 0 and denom != 0:
        return inf.real

    if denom == 0 and denom == 0:
        return nan.real

    return numer/denom


def gain_unilateral(S: tuple[complex, complex, complex, complex], gammaS: int | float | complex, gammaL: int | float | complex) -> float:
    """
    Calculates the unilateral transducer power gain of a 2-port network.

    Args:
        S (tuple[complex, complex, complex, complex]): Unfolded scattering matrix (S11, S21, S12, S22).
        gammaS (int | float | complex): Reflection coefficient of the source.
        gammaL (int | float | complex): Reflection coefficient of the load.

    Returns:
        float: Unilateral transducer power gain of the network in linear scale.
    """

    S11: complex = S[0]
    S21: complex = S[1]
    S22: complex = S[3]

    if S21 == 0:
        return 0

    numer: float = (1-abs(gammaS)**2)*(abs(S21)**2)*(1-abs(gammaL)**2)
    denom: float = abs((1-S11*gammaS)*(1-S22*gammaL))**2

    if denom == 0 and denom != 0:
        return inf.real

    if denom == 0 and denom == 0:
        return nan.real

    return numer/denom


def gain_availabe(S: tuple[complex, complex, complex, complex], gammaS: int | float | complex) -> float:
    """
    Calculates the available power gain of a 2-port network.

    Args:
        S (tuple[complex, complex, complex, complex]): Unfolded scattering matrix (S11, S21, S12, S22).
        gammaS (int | float | complex): Reflection coefficient of the source.

    Returns:
        float: Available power gain of the network in linear scale.
    """

    S11: complex = S[0]
    S21: complex = S[1]

    if S21 == 0:
        return 0

    gammaOut: complex = output_reflection(S, gammaS)

    numer: float = (1-abs(gammaS)**2)*abs(S21)**2
    denom: float = (1-abs(gammaOut)**2)*abs(1-S11*gammaS)**2

    if denom == 0 and denom != 0:
        return inf.real

    if denom == 0 and denom == 0:
        return nan.real

    return numer/denom


def gain_operating(S: tuple[complex, complex, complex, complex], gammaL: int | float | complex) -> float:
    """
    Calculates the operating power gain of a 2-port network.

    Args:
        S (tuple[complex, complex, complex, complex]): Unfolded scattering matrix (S11, S21, S12, S22).
        gammaL (int | float | complex): Reflection coefficient of the load.

    Returns:
        float: Operating power gain of the network in linear scale.
    """

    S21: complex = S[1]
    S22: complex = S[3]

    if S21 == 0:
        return 0

    gammaIn: complex = input_reflection(S, gammaL)

    numer: float = (1-abs(gammaL)**2)*abs(S21)**2
    denom: float = (1-abs(gammaIn)**2)*abs(1-S22*gammaL)**2

    if denom == 0 and denom != 0:
        return inf.real

    if denom == 0 and denom == 0:
        return nan.real

    return numer/denom


def gain_maximum(S: tuple[complex, complex, complex, complex]) -> float:
    """
    Calculates the power gain of a 2-port network when load and source are matched for maximum power transfer. If the K factor is not greater than one, returns NaN.

    Args:
        S (tuple[complex, complex, complex, complex]): Unfolded scattering matrix (S11, S21, S12, S22).

    Returns:
        float: Maximum power gain of the network in linear scale.
    """

    S11: complex = S[0]
    S21: complex = S[1]
    S12: complex = S[2]
    S22: complex = S[3]
    K: float = rollet(S)

    if K <= 1:
        return nan.real

    if S21 == 0:
        return 0

    if S12 == 0:
        if abs(S11) == 1 or abs(S22) == 1:
            return nan.real

        Ggmax: float = abs(1/(1-abs(S11)**2))
        Glmax: float = abs(1/(1-abs(S22)**2))
        return Ggmax*Glmax*abs(S21)**2

    if isinf(S21) and isinf(S12):
        return nan.real

    if isfinite(S21) and S12 == 0:
        return inf.real

    return abs(S21/S12)/(K+mt.sqrt(K**2-1))


def gain_maximum_stable(S: tuple[complex, complex, complex, complex]) -> float:
    """
    Calculates the maximum stable power gain of a conditionaly stable 2-port network.

    Args:
        S (tuple[complex, complex, complex, complex]): Unfolded scattering matrix (S11, S21, S12, S22).

    Returns:
        float: Maximum stable power gain of the network in linear scale.
    """

    S21: complex = S[1]
    S12: complex = S[2]
    K: float = rollet(S)

    if isnan(K):
        return nan.real

    if K >= 1:
        return nan.real

    if S21 == 0:
        return 0

    if S12 == 0 and S21 == 0:
        return nan.real

    if isinf(S12) and isinf(S21):
        return nan.real

    if isinf(S12) and isfinite(S21):
        return inf.real

    return abs(S21/S12)


def unilateral_test(S: tuple[complex, complex, complex, complex]) -> tuple[float, float, float]:
    """
    Calculates the unilaterality factor and inferior and superior limits for the ratio between transducer gain and unilateral transducer gain.

    Args:
        S (tuple[complex, complex, complex, complex]): Unfolded scattering matrix (S11, S21, S12, S22).

    Returns:
        tuple[float, float, float]: (U, lim_inf, lim_sup)
            U (float): Unilaterality factor.
            lim_inf: Inferior limit for G_T/G_TU.
            lim_sup: Superior limit for G_T/G_TU.
    """

    S11: complex = S[0]
    S21: complex = S[1]
    S12: complex = S[2]
    S22: complex = S[3]

    U: float

    numerU: float = abs(S12*S21*S11*S22)
    denomU: float = (1-abs(S11)**2)*(1-abs(S22)**2)

    if (numerU == 0 and denomU == 0) or (isinf(numerU) and isinf(denomU)):
        U = nan.real
    elif numerU != 0 and numerU == 0:
        U = inf.real
    else:
        U = numerU/denomU

    lim_inf: float
    lim_sup: float

    if U == -1:
        lim_inf = inf.real
    else:
        lim_inf = 1/((1+U)**2)

    if U == 1:
        lim_sup = inf.real
    else:
        lim_sup = 1/((1-U)**2)

    return (U, lim_inf, lim_sup)


def x_to_lc(x: int | float, f: int | float) -> float:
    """
    Calculates the equivalent capacitance or inductance from a given reactance. When the choice between capacitance or inductance cannot be made (x=0, x=inf, f=0 or x=inf), returns NaN. If the frequency is negative, the absolute value is used.

    Args:
        x (float): Reactance to be converted.
        f (float): Frequency of operation.

    Returns:
        float: Equivalent capacitance or inductance signified by sign (positive <-> inductance, negative <-> capacitance).
    """

    if x*f == 0 or isinf(x*f):
        return nan.real

    if x > 0:
        return x/(2*pi*abs(f))

    return 1/(2*pi*abs(f)*x)


def l_matching_network(ZS: int | float | complex, ZL: int | float | complex, f: int | float, block_DC: bool = False) -> tuple[float, float, float, bool]:
    """
    Calculate the L matching network that adapts a power source with impedance ZS to a load with impedance ZL at a given frequency. The algorithm adds parallel or series reactances to load and source to compensate their imaginary parts. Then it uses the the simplified resistance transformation equations to design a L matching network. Those added reactances are then incorporated into the network.

    Args:
        ZS (int | float | complex): Source impedance.
        ZL (int | float | complex): Load impedance.
        f (int | float): Frequency of operation.
        block_DC (bool, optional): Optional flag defining if the network is low-pass or high-pass. Defaults to False.

    Returns:
        tuple[float, float, float, bool]: (LC1, LC2, Q, source_parallel)
            LC1 (float): Source side capacitance or inductance signified by sign (positive <-> inductance, negative <-> capacitance).
            LC2 (float): Load side capacitance or inductance signified by sign (positive <-> inductance, negative <-> capacitance).
            Q (float): Unloaded Q factor.
            source_parallel (bool): Flag that indicated if the source side is parrallel.
    """

    GS: float = real(1/ZS)
    BS: float = imag(1/ZS)

    RL: float = real(ZL)
    XL: float = imag(ZL)

    if GS*RL > 1:
        LC1p: float
        LC2p: float
        Qp: float
        LC2p, LC1p, Qp = l_matching_network(ZS=ZL, ZL=ZS, f=f, block_DC=block_DC)[:3]
        return (LC1p, LC2p, Qp, False)

    Q: float = mt.sqrt(1/(GS*RL) - 1)

    B1: float = Q*GS
    X2: float = Q*RL

    if block_DC:
        B1 *= -1
        X2 *= -1

    B1p: float = B1 - BS
    X2p: float = X2 - XL

    LC1: float = x_to_lc(-1/B1p, f)
    LC2: float = x_to_lc(X2p, f)

    return (LC1, LC2, Q, True)


def two_port_match(S: tuple[complex, complex, complex, complex]) -> tuple[complex, complex]:
    """
    Calculates the source and load refletion coefficients that maximize the power gain. If the coefficients are not calculable or have magnitude greater than 1, returns a tuple of NaN.

    Args:
        S (tuple[complex, complex, complex, complex]): Unfolded scattering matrix (S11, S21, S12, S22).

    Returns:
        tuple[complex, complex]: (gammaS, gammaL)
            gammaS (complex): Optimal refelction coefficient on source.
            gammaL (complex): Optimal reflection coefficient on load.
    """

    S11: complex = S[0]
    S21: complex = S[1]
    S12: complex = S[2]
    S22: complex = S[3]
    delta: complex = S11*S22-S12*S21

    B: float = 1 + abs(S11)**2 - abs(S22)**2 - abs(delta)**2
    C: complex = S11 - delta*conj(S22)

    if C == 0 and B == 0:
        return (nan, nan)

    if C == 0 and B != 0:
        return (0, conj(S22))

    gammaS: complex = (B + cm.sqrt(B**2 - 4*abs(C)**2))/(2*C)

    if abs(gammaS) > 1:
        gammaS = (B - cm.sqrt(B**2 - 4*abs(C)**2))/(2*C)

    if abs(gammaS) > 1:
        return (nan, nan)

    gammaL: complex = conj(output_reflection(S, gammaS))

    return (gammaS, gammaL)
