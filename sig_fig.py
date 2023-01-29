import math as mt

si_prefix: dict[int, str] = {-30: " q",
                             -27: " r",
                             -24: " y",
                             -21: " z",
                             -18: " a",
                             -15: " f",
                             -12: " p",
                             -9: " n",
                             -6: " μ",
                             -3: " m",
                             0: " ",
                             3: " k",
                             6: " M",
                             9: " G",
                             12: " T",
                             15: " P",
                             18: " E",
                             21: " Z",
                             24: " Y",
                             27: " R",
                             30: " Q"}


def round_sci(x: int | float, n: int = 5, unit: str = "") -> str:
    """
    Function that rounds and format a value in scientific notation.

    Args:
        x (int | float): Value to be formatted.
        n (int, optional): Amount of significant figures. Defaults to 5.
        unit (str, optional): Unit of the value. Defaults to "".

    Returns:
        str: Formatted string.
    """

    if unit != "":
        unit = " "+unit

    if n == 0:
        return ""

    if x < 0:
        return "-"+round_sci(x=abs(x), n=n, unit=unit)

    if mt.isinf(x):
        return "inf"+unit

    if mt.isnan(x):
        return "NaN"

    if x == 0:
        long_string = "0."+(n-1)*"0"
        if long_string[-1] == ".":
            long_string = long_string[:-1]
        return long_string + unit

    sci_exp = mt.floor(mt.log10(x))
    x_norm = x/(10**sci_exp)
    whole, frac = str(x_norm).split(".")
    long_string = whole + frac

    if len(long_string) == n:
        if sci_exp == 0:
            return whole + "." + frac + unit
        elif sci_exp == 1:
            return whole + "." + frac + "×10" + unit
        else:
            return whole + "." + frac + "×10^" + str(sci_exp) + unit
    elif len(long_string) < n:
        long_string += (n - len(long_string) + 1) * "0"

    long_list = [int(char) for char in long_string]

    digit_before_trunc = long_list[n-1]
    digit_after_trunc = long_list[n]

    if digit_after_trunc > 5:
        long_list[n-1] += 1
    elif digit_after_trunc == 5:
        if set(long_list[n:]).intersection(set([1, 2, 3, 4, 5, 6, 7, 8, 9])):
            long_list[n-1] += 1
        elif digit_before_trunc % 2 == 1:
            long_list[n-1] += 1

    long_list = long_list[:n]

    index = len(long_list) - 1
    while index > 0:
        current_digit = long_list[index]
        if current_digit > 9:
            long_list[index] -= 10
            long_list[index-1] += 1
        index -= 1

    if long_list[0] > 9:
        sci_exp += 1
        long_list[0] -= 10
        long_list = [1] + long_list

    long_list = long_list[:n]

    long_string = "".join([str(digit) for digit in long_list])
    whole = long_string[:1]
    frac = long_string[1:]
    if sci_exp == 0:
        return whole + "." + frac + unit
    elif sci_exp == 1:
        return whole + "." + frac + "×10" + unit
    else:
        return whole + "." + frac + "×10^" + str(sci_exp) + unit


def round_fix(x: int | float, n: int = 5, unit: str = "") -> str:
    """
    Function that rounds and format a value in standard fixed point notation.

    Args:
        x (int | float): Value to be formatted.
        n (int, optional): Amount of significant figures. Defaults to 5.
        unit (str, optional): Unit of the value. Defaults to "".

    Returns:
        str: Formatted string.
    """

    if unit != "":
        unit = " "+unit

    if n == 0:
        return ""

    if mt.isnan(x):
        return "NaN"

    if x < 0:
        return "-"+round_fix(x=abs(x), n=n, unit=unit)

    if mt.isinf(x):
        return "inf"+unit

    sci_str = round_sci(x=x, n=n, unit="")

    if "×" not in sci_str:
        return sci_str + unit

    if "^" not in sci_str:
        num_exp = 1
    else:
        num_exp = int(sci_str.split("^")[1])

    long_str = sci_str.split("×")[0]

    if "." not in long_str:
        long_str += "."

    while num_exp != 0:
        dot_index = long_str.index(".")
        num_len = len(long_str)

        if num_exp > 0:
            if dot_index == num_len-1:
                long_str = long_str[:-1] + "0" + "."
            else:
                char_after_dot = long_str[dot_index+1]
                long_str = long_str[:dot_index] + char_after_dot + "." + long_str[dot_index+2:]
            num_exp -= 1
        else:
            if dot_index == 0:
                long_str = ".0" + long_str[1:]
            elif dot_index == 1:
                char_before_dot = long_str[0]
                long_str = "0." + char_before_dot + long_str[2:]
            else:
                char_before_dot = long_str[dot_index-1]
                long_str = long_str[:dot_index-2] + "." + char_before_dot + long_str
            num_exp += 1

    if long_str[-1] == ".":
        long_str = long_str[:-1]

    return long_str + unit


def round_eng(x: int | float, n: int = 5, unit: str = "") -> str:
    """
    Function that rounds and format a value in enginerring notation (using the appropriate SI prefix).

    Args:
        x (int | float): Value to be formatted.
        n (int, optional): Amount of significant figures. Defaults to 5.
        unit (str, optional): Unit of the value. Defaults to "".

    Returns:
        str: Formatted string.
    """

    if n == 0:
        return ""

    if x < 0:
        return "-"+round_eng(x=abs(x), n=n, unit=unit)

    if mt.isinf(x):
        if unit == "":
            return "inf"
        else:
            return "inf "+unit

    if mt.isnan(x):
        return "NaN"

    if x == 0:
        long_string = "0."+(n-1)*"0"
        if long_string[-1] == ".":
            long_string = long_string[:-1]
        if unit == "":
            return long_string
        else:
            return long_string + " " + unit

    si_exp = 3*mt.floor(mt.log10(x)/3)
    x_norm = x/(10**si_exp)
    whole, frac = str(x_norm).split(".")
    decimal_point = len(whole)
    long_string = whole + frac

    if len(long_string) == n:
        return whole + "." + frac + si_prefix[si_exp] + unit
    elif len(long_string) < n:
        long_string += (n - len(long_string) + 1) * "0"

    long_list = [int(char) for char in long_string]

    digit_before_trunc = long_list[n-1]
    digit_after_trunc = long_list[n]

    if digit_after_trunc > 5:
        long_list[n-1] += 1
    elif digit_after_trunc == 5:
        if set(long_list[n:]).intersection(set([1, 2, 3, 4, 5, 6, 7, 8, 9])):
            long_list[n-1] += 1
        elif digit_before_trunc % 2 == 1:
            long_list[n-1] += 1

    long_list = long_list[:n]

    index = len(long_list) - 1
    while index > 0:
        current_digit = long_list[index]
        if current_digit > 9:
            long_list[index] -= 10
            long_list[index-1] += 1
        index -= 1

    if long_list[0] > 9:
        decimal_point += 1
        long_list[0] -= 10
        long_list = [1] + long_list

    if decimal_point > 3:
        si_exp += 3
        decimal_point -= 3

    long_string = "".join([str(digit) for digit in long_list])

    n_string = len(long_string)
    if n_string < decimal_point:
        long_string += (decimal_point - n_string) * "0"

    long_string = long_string[:decimal_point] + "." + long_string[decimal_point:]

    whole, frac = long_string.split(".")
    if len(whole) >= n:
        long_string = whole
    else:
        n_rest = n - len(whole)
        long_string = whole + "." + frac[:n_rest]

    return long_string + si_prefix[si_exp] + unit


def complex_round_sci(z: int | float | complex, n: int = 5, unit: str = "") -> str:
    """
    Function that rounds and format a value as a complex number with real and imaginary parts in scientific notation.

    Args:
        z (int | float | complex): Value to be formatted.
        n (int, optional): Amount of significant figures. Defaults to 5.
        unit (str, optional): Unit of the value. Defaults to "".

    Returns:
        str: Formatted string.
    """

    if unit != "":
        unit = " "+unit

    a = complex(z).real
    b = complex(z).imag

    if mt.isinf(a) or mt.isinf(b):
        return "inf"+unit

    if a == 0 and b == 0:
        return round_sci(x=0, n=n, unit=unit)
    if b == 0:
        return round_sci(x=a, n=n, unit=unit)

    if a == 0 and b > 0:
        return "j"+round_sci(abs(b), n=n)
    elif a == 0 and b < 0:
        return "-j"+round_sci(abs(b), n=n)

    if b < 0:
        b_str = "-j" + round_sci(abs(b), n=n)
    else:
        b_str = "+j" + round_sci(abs(b), n=n)

    a_str = round_sci(x=a, n=n)

    if unit == "":
        return a_str + b_str

    return "(" + a_str + b_str + ")" + unit


def complex_round_fix(z: int | float | complex, n: int = 5, unit: str = "") -> str:
    """
    Function that rounds and format a value as a complex number with real and imaginary parts in standard fixed point notation.

    Args:
        z (int | float | complex): Value to be formatted.
        n (int, optional): Amount of significant figures. Defaults to 5.
        unit (str, optional): Unit of the value. Defaults to "".

    Returns:
        str: Formatted string.
    """

    if unit != "":
        unit = " "+unit

    a = complex(z).real
    b = complex(z).imag

    if mt.isinf(a) or mt.isinf(b):
        return "inf"+unit

    if a == 0 and b == 0:
        return round_fix(x=0, n=n, unit=unit)
    if b == 0:
        return round_fix(x=a, n=n, unit=unit)

    if a == 0 and b > 0:
        return "j"+round_fix(abs(b), n=n)
    elif a == 0 and b < 0:
        return "-j"+round_fix(abs(b), n=n)

    if b < 0:
        b_str = "-j" + round_fix(abs(b), n=n)
    else:
        b_str = "+j" + round_fix(abs(b), n=n)

    a_str = round_fix(x=a, n=n)

    if unit == "":
        return a_str + b_str

    return "(" + a_str + b_str + ")" + unit
