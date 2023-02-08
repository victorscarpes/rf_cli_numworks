import cmath as cm
import sig_fig as sf
from rf_active_tools import *

print("Enter operating frequency in MHz.")

str_in: str = input("f: ")
f: float = eval(str_in)*1e6

print("\nEnter S parameters as linear magnitude and")
print("phase in degrees separated by comma with no")
print("spaces. If a single value is entered")
print("without any comma, it is considered as the")
print("magnitude and the phase defaults to zero.")

str_in = input("\nS_11: ")
s11_mag_str: str = str_in
s11_deg_str: str = "0"
if "," in str_in:
    s11_mag_str = str_in.split(",")[0]
    s11_deg_str = str_in.split(",")[1]
s11_mag: float = eval(s11_mag_str)
s11_rad: float = eval(s11_deg_str)*pi/180
S11: complex = pol(s11_mag, s11_rad)
del s11_mag, s11_mag_str, s11_rad, s11_deg_str

str_in = input("S_21: ")
s21_mag_str: str = str_in
s21_deg_str: str = "0"
if "," in str_in:
    s21_mag_str = str_in.split(",")[0]
    s21_deg_str = str_in.split(",")[1]
s21_mag: float = eval(s21_mag_str)
s21_rad: float = eval(s21_deg_str)*pi/180
S21: complex = pol(s21_mag, s21_rad)
del s21_mag, s21_mag_str, s21_rad, s21_deg_str

str_in = input("S_12: ")
s12_mag_str: str = str_in
s12_deg_str: str = "0"
if "," in str_in:
    s12_mag_str = str_in.split(",")[0]
    s12_deg_str = str_in.split(",")[1]
s12_mag: float = eval(s12_mag_str)
s12_rad: float = eval(s12_deg_str)*pi/180
S12: complex = pol(s12_mag, s12_rad)
del s12_mag, s12_mag_str, s12_rad, s12_deg_str

str_in = input("S_22: ")
s22_mag_str: str = str_in
s22_deg_str: str = "0"
if "," in str_in:
    s22_mag_str = str_in.split(",")[0]
    s22_deg_str = str_in.split(",")[1]
s22_mag: float = eval(s22_mag_str)
s22_rad: float = eval(s22_deg_str)*pi/180
S22: complex = pol(s22_mag, s22_rad)
del s22_mag, s22_mag_str, s22_rad, s22_deg_str

S = (S11, S21, S12, S22)


print("\nEnter reference, load and source impedances")
print("in ohms. Values default to 50 Ω when not")
print("specified.")

str_in = input("\nZ_0: ")
Z0: complex = complex(50)
if str_in != "":
    Z0 = eval(str_in)

str_in = input("Z_S: ")
ZS: complex = complex(50)
if str_in != "":
    ZS = eval(str_in)

str_in = input("Z_L: ")
ZL: complex = complex(50)
if str_in != "":
    ZL = eval(str_in)

gammaS: complex = z_to_s(ZS, Z0)
gammaL: complex = z_to_s(ZL, Z0)

Gt: float = gain_transducer(S, gammaS, gammaL)
Ga: float = gain_availabe(S, gammaS)
Go: float = gain_operating(S, gammaL)
Gmax: float = gain_maximum(S)
Gstab: float = gain_maximum_stable(S)
Gtu: float = gain_unilateral(S, gammaS, gammaL)

K: float = rollet(S)
delta: complex = S11*S22 - S12*S21
muL: float
muS: float
muS, muL = mu_stab(S)

Ocs: complex
rs: float
Ocs, rs = source_stab_circle(S)

Ocl: complex
rl: float
Ocl, rl = load_stab_circle(S)

U: float
lim_inf: float
lim_sup: float
U, lim_inf, lim_sup = unilateral_test(S)
gain_error: float = 100*max(abs(1-lim_inf), abs(1-lim_sup))

gammaS_match: complex
gammaL_match: complex
gammaS_match, gammaL_match = two_port_match(S)
ZS_match: complex = s_to_z(gammaS_match, Z0)
ZL_match: complex = s_to_z(gammaL_match, Z0)

gamma_in: complex = input_reflection(S, gammaL)
gamma_out: complex = output_reflection(S, gammaS)
Zin: complex = s_to_z(gamma_in, Z0)
Zout: complex = s_to_z(gamma_out, Z0)

Zin_opt: complex = conj(ZS_match)
Zout_opt: complex = conj(ZL_match)

LC1_in_lp: float
LC2_in_lp: float
Q_in_lp: float
source_parallel_lp: bool
LC1_in_lp, LC2_in_lp, Q_in_lp, source_parallel_lp = l_matching_network(ZS, Zin_opt, f)

LC1_in_hp: float
LC2_in_hp: float
Q_in_hp: float
source_parallel_hp: bool
LC1_in_hp, LC2_in_hp, Q_in_hp, source_parallel_hp = l_matching_network(ZS, Zin_opt, f, True)

LC1_out_lp: float
LC2_out_lp: float
Q_out_lp: float
amp_parallel_lp: bool
LC1_out_lp, LC2_out_lp, Q_out_lp, amp_parallel_lp = l_matching_network(Zout_opt, ZL, f)

LC1_out_hp: float
LC2_out_hp: float
Q_out_hp: float
amp_parallel_hp: bool
LC1_out_hp, LC2_out_hp, Q_out_hp, amp_parallel_hp = l_matching_network(Zout_opt, ZL, f, True)

loop_flag: bool = True


while loop_flag:
    print("\n"+43*"-")
    print("0 - Gain metrics")
    print("1 - Unilaterality metrics")
    print("2 - Stability metrics")
    print("3 - Stability circles")
    print("4 - Simultaneous matching")
    print("5 - Input low-pass L network")
    print("6 - Input high-pass L network")
    print("7 - Output low-pass L network")
    print("8 - Output high-pass L network")
    print("9 - Input and output impedances")

    print("\nEnter operation or press [ENTER] to end")
    str_in = input("the program: ")

    if str_in == "0":  # Gain metrics
        print("\n"+43*"-")
        print("G_T = "+sf.round_fix(dB(Gt), unit="dB"))
        print("G_A = "+sf.round_fix(dB(Ga), unit="dB"))
        print("G_OP = "+sf.round_fix(dB(Go), unit="dB"))
        print("G_TU = "+sf.round_fix(dB(Gtu), unit="dB"))

        if K > 1 and abs(delta) < 1:
            print("G_max = "+sf.round_fix(dB(Gmax), unit="dB"))
        elif abs(K) < 1:
            print("G_stab = "+sf.round_fix(dB(Gstab), unit="dB"))

        input("\nPress [ENTER] to continue.")

    elif str_in == "1":  # Unilaterality metrics
        print("\n"+43*"-")
        print("G_TU = "+sf.round_fix(dB(Gtu), unit="dB"))
        print("U = "+sf.round_fix(U))
        print(sf.round_fix(lim_inf)+" < G_T/G_TU < "+sf.round_fix(lim_sup))
        print("Δ% = "+sf.round_fix(gain_error)+"%")

        input("\nPress [ENTER] to continue.")

    elif str_in == "2":  # Stability metrics
        print("\n"+43*"-")
        print("K = "+sf.round_fix(K))
        print("|Δ| = "+sf.round_fix(abs(delta)))
        print("arg(Δ) = "+sf.round_fix(cm.phase(delta)*180/pi, unit="deg"))
        print("μ_S = "+sf.round_fix(muS))
        print("μ_L = "+sf.round_fix(muL))

        if abs(K) == 1:
            print("\nNetwork is unmatchable.")
        elif K > 1 and abs(delta) < 1:
            print("\nNetwork is matchable and unconditionally")
            print("stable.")
        elif K > 1 and abs(delta) > 1:
            print("\nNetwork is matchable and conditionally")
            print("stable.")
        elif K < -1:
            print("\nNetwork is unmatchable and unconditionally")
            print("unstable.")

        input("\nPress [ENTER] to continue.")

    elif str_in == "3":  # Stability circles
        print("\n"+43*"-")

        print("|O_S| = "+sf.round_fix(abs(Ocs)))
        print("arg(O_S) = "+sf.round_fix(cm.phase(Ocs)*180/pi, unit="deg"))
        print("r_S = "+sf.round_fix(rs))

        print("\n|O_L| = "+sf.round_fix(abs(Ocl)))
        print("arg(O_L) = "+sf.round_fix(cm.phase(Ocl)*180/pi, unit="deg"))
        print("r_L = "+sf.round_fix(rl))

        input("\nPress [ENTER] to continue.")

    elif str_in == "4":  # Simultaneous matching
        print("\n"+43*"-")

        if isnan(gammaS_match) or isnan(gammaL_match):
            print("Unable to match network.")

        else:
            print("|Γ_S| = "+sf.round_fix(abs(gammaS_match)))
            print("arg(Γ_S) = "+sf.round_fix(cm.phase(gammaS_match)*180/pi, unit="deg"))

            print("\n|Γ_L| = "+sf.round_fix(abs(gammaL_match)))
            print("arg(Γ_L) = "+sf.round_fix(cm.phase(gammaL_match)*180/pi, unit="deg"))

            print("\nZ_S = "+sf.complex_round_fix(ZS_match, unit="Ω"))
            print("Z_L = "+sf.complex_round_fix(ZL_match, unit="Ω"))

        input("\nPress [ENTER] to continue.")

    elif str_in == "5":  # Input low-pass L network
        print("\n"+43*"-")

        print("Q = "+sf.round_fix(Q_in_lp))

        if source_parallel_lp:
            print("\nIn parallel with source:")
        else:
            print("\nIn series with source:")

        if LC1_in_lp > 0:
            print("L = "+sf.round_eng(LC1_in_lp, unit="H"))
        else:
            print("C = "+sf.round_eng(-LC1_in_lp, unit="F"))

        if source_parallel_lp:
            print("\nIn series with amp input:")
        else:
            print("\nIn parallel with amp input:")

        if LC2_in_lp > 0:
            print("L = "+sf.round_eng(LC2_in_lp, unit="H"))
        else:
            print("C = "+sf.round_eng(-LC2_in_lp, unit="F"))

        input("\nPress [ENTER] to continue.")

    elif str_in == "6":  # Input high-pass L network
        print("\n"+43*"-")

        print("Q = "+sf.round_fix(Q_in_hp))

        if source_parallel_hp:
            print("\nIn parallel with source:")
        else:
            print("\nIn series with source:")

        if LC1_in_hp > 0:
            print("L = "+sf.round_eng(LC1_in_hp, unit="H"))
        else:
            print("C = "+sf.round_eng(-LC1_in_hp, unit="F"))

        if source_parallel_hp:
            print("\nIn series with amp input:")
        else:
            print("\nIn parallel with amp input:")

        if LC2_in_hp > 0:
            print("L = "+sf.round_eng(LC2_in_hp, unit="H"))
        else:
            print("C = "+sf.round_eng(-LC2_in_hp, unit="F"))

        input("\nPress [ENTER] to continue.")

    elif str_in == "7":  # Output low-pass L network
        print("\n"+43*"-")

        print("Q = "+sf.round_fix(Q_out_lp))

        if amp_parallel_lp:
            print("\nIn parallel with amp output:")
        else:
            print("\nIn series with amp output:")

        if LC1_out_lp > 0:
            print("L = "+sf.round_eng(LC1_out_lp, unit="H"))
        else:
            print("C = "+sf.round_eng(-LC1_out_lp, unit="F"))

        if amp_parallel_lp:
            print("\nIn series with load:")
        else:
            print("\nIn parallel with load:")

        if LC2_out_lp > 0:
            print("L = "+sf.round_eng(LC2_out_lp, unit="H"))
        else:
            print("C = "+sf.round_eng(-LC2_out_lp, unit="F"))

        input("\nPress [ENTER] to continue.")

    elif str_in == "8":  # Output high-pass L network
        print("\n"+43*"-")

        print("Q = "+sf.round_fix(Q_out_hp))

        if amp_parallel_hp:
            print("\nIn parallel with amp output:")
        else:
            print("\nIn series with amp output:")

        if LC1_out_hp > 0:
            print("L = "+sf.round_eng(LC1_out_hp, unit="H"))
        else:
            print("C = "+sf.round_eng(-LC1_out_hp, unit="F"))

        if amp_parallel_hp:
            print("\nIn series with load:")
        else:
            print("\nIn parallel with load:")

        if LC2_out_hp > 0:
            print("L = "+sf.round_eng(LC2_out_hp, unit="H"))
        else:
            print("C = "+sf.round_eng(-LC2_out_hp, unit="F"))

        input("\nPress [ENTER] to continue.")

    elif str_in == "9":  # Input and output impedances
        print("\n"+43*"-")

        print("|Γ_in| = "+sf.round_fix(abs(gamma_in)))
        print("arg(Γ_in) = "+sf.round_fix(cm.phase(gamma_in)*180/pi, unit="deg"))

        print("\n|Γ_out| = "+sf.round_fix(abs(gamma_out)))
        print("arg(Γ_out) = "+sf.round_fix(cm.phase(gamma_out)*180/pi, unit="deg"))

        print("\nZ_in = "+sf.complex_round_fix(Zin, unit="Ω"))
        print("Z_out = "+sf.complex_round_fix(Zout, unit="Ω"))

        input("\nPress [ENTER] to continue.")

    else:  # End program
        loop_flag = False
