"""WBGT calculator for Robert Street Park (Team 207).

WBGT = 0.7Tnw + 0.2Tg + 0.1Ta.
Run: python3 wbgt.py. Sources are listed in README.md.
"""
import math
import random

SIG = 5.67e-8          # Stefan-Boltzmann constant, W/m²K⁴
K = 273.15

# Design condition
TA = 31.0              # Heat Warning air temperature, °C
HUMIDEX = 40.0         # Heat Warning humidex
LAT = 43.66            # Toronto latitude, degrees
DECL = 21.5            # Mid-July solar declination, degrees
HOUR = 15.0            # Hour angle, about 2:20 pm EDT
F_DIR, F_DIF = 0.67, 0.33            # Direct / diffuse sunlight share
U10 = 12.9 / 3.6       # July mean wind at 10 m, m/s
Z0 = 0.5               # Roughness length for a treed park, m
Z_GLOBE = 1.1          # Globe height for a standing person, m
ALB_GROUND = 0.20      # Grass albedo
ALB_GLOBE, EPS_GLOBE = 0.05, 0.95    # Black globe properties
EPS_GROUND = 0.95

# Team readings, Oct 12 2026 (CDS Table A2)
GROUND_SUN_OCT = (13.5 + 5.6 + 4.2 + 9.3) / 4      # Ground above air in sun
GROUND_SHADE_OCT = (9.8 + 2.3 - 0.2 + 0.2) / 4     # Ground above air in shade

# Design inputs
TAU_FABRIC = 0.10      # AD01 membrane solar transmission
ABS_FABRIC = 0.22      # AD01 membrane solar absorption
ABS_PANEL = 0.35       # AD02/AD03 light panel solar absorption
ETA_MIST = 0.20        # AD02 mist: share of possible cooling
DT_TUBES = 2.0         # AD03 terracotta tube air cooling, °C


# Psychrometrics and sun
def vapour_pressure_from_humidex(T, H):
    """Vapour pressure (hPa) and dew point from humidex."""
    e = (H - T) / 0.5555 + 10
    g = math.log(e / 6.112)
    return e, 237.7 * g / (17.27 - g)


def rh(T, e):
    return 100 * e / (6.112 * math.exp(17.67 * T / (T + 243.5)))


def wetbulb(T, RH):
    """Stull (2011) wet-bulb temperature, °C."""
    return (T * math.atan(0.151977 * (RH + 8.313659) ** 0.5) + math.atan(T + RH) - math.atan(RH - 1.676331)
            + 0.00391838 * RH ** 1.5 * math.atan(0.023101 * RH) - 4.686035)


def cos_zenith(decl, hour):
    p, d, h = map(math.radians, (LAT, decl, hour))
    return math.sin(p) * math.sin(d) + math.cos(p) * math.cos(d) * math.cos(h)


def haurwitz(cosz):
    """Haurwitz (1945) clear-sky sunlight, W/m²."""
    return 1098 * cosz * math.exp(-0.057 / cosz) if cosz > 0 else 0.0


def wind_at(z):
    """Log wind profile from the 10 m wind."""
    return U10 * math.log(z / Z0) / math.log(10 / Z0)


E_A, TD = vapour_pressure_from_humidex(TA, HUMIDEX)
RH0 = rh(TA, E_A)
TW = wetbulb(TA, RH0)
Z = math.acos(cos_zenith(DECL, HOUR))
S_GLOBAL = haurwitz(math.cos(Z))                    # July design sunlight
OCT_SUN = haurwitz(cos_zenith(-7.5, 12.0))          # Sunlight when the team measured
D_SUN = GROUND_SUN_OCT * S_GLOBAL / OCT_SUN         # Sunlit ground excess, scaled to July
D_SHADE = GROUND_SHADE_OCT
U = wind_at(Z_GLOBE)
EPS_SKY = 0.575 * E_A ** (1 / 7)                    # Sky emissivity (Hunter and Minyard)
H_C = 0.115 * (U * 3600) ** 0.58                    # Globe convection coefficient, W/m²K


# Geometry
def sky_view(half_w, length, height, n=200000, seed=1):
    """Open-sky share under a rectangular canopy (Monte Carlo)."""
    rnd = random.Random(seed)
    open_ = 0
    for _ in range(n):
        cz = rnd.random()
        phi = rnd.random() * 2 * math.pi
        sz = math.sqrt(1 - cz * cz)
        x, y = sz * math.cos(phi) / cz * height, sz * math.sin(phi) / cz * height
        if not (abs(x) <= half_w and abs(y) <= length / 2):
            open_ += 1
    return open_ / n


FSKY_AD01 = sky_view(1.5, 200, 3.0 - Z_GLOBE)     # 3 m wide ribbon, 3 m high
FSKY_AD02 = sky_view(3.0, 6.0, 3.0 - Z_GLOBE)     # 6 m x 6 m node, 3 m high
FSKY_AD03 = sky_view(2.0, 5.0, 2.7 - Z_GLOBE)     # 5 m x 4 m module, 2.7 m high


# Heat balances
def canopy_temp(absorb, Ta=TA):
    """Sunlit canopy temperature from its heat balance."""
    h = 5.7 + 3.8 * wind_at(3.0)
    L_sky = EPS_SKY * SIG * (Ta + K) ** 4
    L_gr = EPS_GROUND * SIG * (Ta + D_SHADE + K) ** 4
    lo, hi = Ta, Ta + 40
    for _ in range(60):
        T = (lo + hi) / 2
        loss = 2 * h * (T - Ta) + 0.9 * (2 * SIG * (T + K) ** 4 - L_sky - L_gr)
        lo, hi = (T, hi) if loss < absorb * S_GLOBAL else (lo, T)
    return T


T_FABRIC = canopy_temp(ABS_FABRIC)
T_PANEL = canopy_temp(ABS_PANEL)


def globe(Ta, tau_dir, fsky, tau_canopy, t_canopy, t_ground):
    """Globe temperature Tg (°C) from its heat balance."""
    Sdir, Sdif = S_GLOBAL * F_DIR, S_GLOBAL * F_DIF
    sw_dir = tau_dir * Sdir / (4 * math.cos(Z))
    sw_dif = 0.5 * Sdif * (fsky + (1 - fsky) * tau_canopy)
    s_local = tau_dir * Sdir + Sdif * (fsky + (1 - fsky) * tau_canopy)     # Sunlight reaching this spot
    sw_ref = 0.5 * ALB_GROUND * s_local
    Q_sw = (1 - ALB_GLOBE) * (sw_dir + sw_dif + sw_ref)
    L_sky = EPS_SKY * SIG * (Ta + K) ** 4
    L_can = 0.90 * SIG * (t_canopy + K) ** 4
    L_gr = EPS_GROUND * SIG * (t_ground + K) ** 4
    Q_lw = EPS_GLOBE * (0.5 * (fsky * L_sky + (1 - fsky) * L_can) + 0.5 * L_gr)
    lo, hi = Ta - 10, Ta + 60
    for _ in range(80):
        Tg = (lo + hi) / 2
        out = EPS_GLOBE * SIG * (Tg + K) ** 4 + H_C * (Tg - Ta)
        lo, hi = (Tg, hi) if out < Q_sw + Q_lw else (lo, Tg)
    return Tg, dict(sw_dir=sw_dir, sw_dif=sw_dif, sw_ref=sw_ref, Q_sw=Q_sw, Q_lw=Q_lw, S_local=s_local)


def natural_wb(Tw, S_local):
    """Hunter and Minyard Eq. 3, solar coefficient 0.0021."""
    return min(Tw + 0.0021 * S_local - 0.42 * U + 1.93, TA)


def case(name, tau_dir=1.0, fsky=1.0, tau_canopy=1.0, t_canopy=None, ground_excess=None, dTa_evap=0.0):
    """WBGT at one location; evaporation keeps wet-bulb unchanged."""
    Ta = TA - dTa_evap
    e = E_A
    if dTa_evap:
        lo, hi = E_A, 60.0
        for _ in range(60):
            e = (lo + hi) / 2
            lo, hi = (e, hi) if wetbulb(Ta, rh(Ta, e)) < TW else (lo, e)
    ge = D_SUN if ground_excess is None else ground_excess
    Tg, parts = globe(Ta, tau_dir, fsky, tau_canopy, Ta if t_canopy is None else t_canopy, Ta + ge)
    Tn = natural_wb(TW, parts['S_local'])
    return dict(name=name, Ta=Ta, RH=rh(Ta, e), Tw=TW, S=parts['S_local'], Tn=Tn, Tg=Tg,
                WBGT=0.7 * Tn + 0.2 * Tg + 0.1 * Ta, parts=parts)


DT_MIST = ETA_MIST * (TA - TW)


def ad01():
    return case('AD01 under the Cooling Ribbon', tau_dir=TAU_FABRIC, fsky=FSKY_AD01, tau_canopy=TAU_FABRIC,
                t_canopy=T_FABRIC, ground_excess=D_SHADE)


CASES = [
    case('Open sun (existing loop)'),
    case('Existing tree shade', tau_dir=0.0, fsky=0.3, tau_canopy=0.05, t_canopy=TA + 1, ground_excess=D_SHADE),
    ad01(),
    case('AD02 under a node canopy, mist on', tau_dir=0.0, fsky=FSKY_AD02, tau_canopy=0.0, t_canopy=T_PANEL,
         ground_excess=D_SHADE, dTa_evap=DT_MIST),
    case('AD02 on the loop, mist only', dTa_evap=DT_MIST),
    case('AD03 inside a module', tau_dir=0.0, fsky=FSKY_AD03, tau_canopy=0.0, t_canopy=T_PANEL,
         ground_excess=D_SHADE, dTa_evap=DT_TUBES),
]

# Effect of each design feature, step by step
STEPS = {
    'AD01': [('Open sun', CASES[0]), ('+ fabric canopy (10% of sun passes)', CASES[2])],
    'AD02': [('Open sun', CASES[0]),
             ('+ opaque node canopy, mist off', case('n', tau_dir=0.0, fsky=FSKY_AD02, tau_canopy=0.0, t_canopy=T_PANEL,
                                                     ground_excess=D_SHADE)),
             ('+ mist ring on', CASES[3])],
    'AD02L': [('Open sun', CASES[0]), ('+ mist on the loop, no canopy', CASES[4])],
    'AD03': [('Open sun', CASES[0]),
             ('+ opaque module canopy, tubes dry', case('m', tau_dir=0.0, fsky=FSKY_AD03, tau_canopy=0.0, t_canopy=T_PANEL,
                                                        ground_excess=D_SHADE)),
             ('+ wet tubes (air 2°C cooler)', CASES[5])],
}


def run_sens():
    """WBGT of open sun and AD01, one input changed."""
    global S_GLOBAL, U, H_C, TAU_FABRIC
    base = (S_GLOBAL, U, TAU_FABRIC)
    out = []
    for label, s, u, tau in [('Base case', *base), ('Hazy sun, 80% of clear sky', 0.8 * base[0], base[1], base[2]),
                             ('Sun 1,000 W/m²', 1000.0, base[1], base[2]), ('Wind 0.5 m/s', base[0], 0.5, base[2]),
                             ('Wind 2 m/s', base[0], 2.0, base[2]), ('Fabric passes 20% of sun', base[0], base[1], 0.20)]:
        S_GLOBAL, U, TAU_FABRIC = s, u, tau
        H_C = 0.115 * (U * 3600) ** 0.58
        out.append((label, case('sun')['WBGT'], ad01()['WBGT']))
    S_GLOBAL, U, TAU_FABRIC = base
    H_C = 0.115 * (U * 3600) ** 0.58
    return out


SENS = run_sens()


# OB02: share of user area below 31 °C
def shaded_share(length, width, shift, n=3600):
    """Average area still shaded when the shadow shifts."""
    tot = 0.0
    for i in range(n):
        a = 2 * math.pi * i / n
        dx, dy = abs(shift * math.cos(a)), abs(shift * math.sin(a))
        tot += max(length - dx, 0) * max(width - dy, 0) / (length * width)
    return tot / n


LOOP_L, LOOP_W = 155.0, 3.0           # Walking loop length and band width, m
PLAY_L, PLAY_W = 16.0, 11.0           # Playground size, m
USER_AREA = LOOP_L * LOOP_W + PLAY_L * PLAY_W
OFFSET_WORST = (3.0 - Z_GLOBE) * math.tan(Z)          # Ribbon shadow shift, sun across path
OFFSET_MEAN = OFFSET_WORST * 2 / math.pi              # Averaged over loop directions
LOOP_SHADED = (LOOP_W - OFFSET_MEAN) / LOOP_W
PLAY_SHADED = shaded_share(PLAY_L, PLAY_W, (4.9 - Z_GLOBE) * math.tan(Z))   # Playground canopy on 4.9 m posts
AD01_COVER = (LOOP_L * LOOP_W * LOOP_SHADED + PLAY_L * PLAY_W * PLAY_SHADED) / USER_AREA
AD02_COVER = 5 * 36.0 / USER_AREA                     # Five 6 m x 6 m node canopies
AD03_COVER = 3 * 20.0 / (USER_AREA + 60.0)            # Three 5 m x 4 m modules


# Report
if __name__ == '__main__':
    print('WBGT calculator - Robert Street Park, Heat Warning threshold\n')
    print(f'Air {TA:.1f} °C, humidex {HUMIDEX:.0f}: e = {E_A:.1f} hPa, dew point {TD:.1f} °C, RH {RH0:.0f}%, wet-bulb {TW:.1f} °C')
    print(f'Sun: zenith {math.degrees(Z):.1f}°, clear-sky {S_GLOBAL:.0f} W/m² (Oct 12 measurement time: {OCT_SUN:.0f} W/m²)')
    print(f'Wind: {U10:.1f} m/s at 10 m -> {U:.2f} m/s at {Z_GLOBE} m; globe convection {H_C:.1f} W/m²K; sky emissivity {EPS_SKY:.2f}')
    print(f'Ground above air: sun +{D_SUN:.1f} °C, shade +{D_SHADE:.1f} °C; canopy above air: fabric +{T_FABRIC - TA:.1f} °C, '
          f'panels +{T_PANEL - TA:.1f} °C')
    print(f'Open sky under canopy: AD01 {FSKY_AD01:.2f}, AD02 {FSKY_AD02:.2f}, AD03 {FSKY_AD03:.2f}; mist cooling {DT_MIST:.1f} °C\n')
    print(f"{'Location':36s} {'Ta':>5s} {'RH%':>4s} {'Tnw':>5s} {'Tg':>5s} {'WBGT':>5s}  CN01 (< 31 °C)")
    for c in CASES:
        print(f"{c['name']:36s} {c['Ta']:5.1f} {c['RH']:4.0f} {c['Tn']:5.1f} {c['Tg']:5.1f} {c['WBGT']:5.1f}  "
              f"{'meets' if c['WBGT'] < 31 else 'fails'}")
    print('\nStep by step (change in WBGT = 0.7 dTnw + 0.2 dTg + 0.1 dTa)')
    for k, st in STEPS.items():
        for (_, p), (lab, c) in zip(st, st[1:]):
            print(f"  {k:5s} {lab:38s} 0.7x{c['Tn'] - p['Tn']:+.2f} + 0.2x{c['Tg'] - p['Tg']:+.1f} + 0.1x{c['Ta'] - p['Ta']:+.1f}"
                  f" = {c['WBGT'] - p['WBGT']:+.2f}  -> {c['WBGT']:.1f} °C")
    print(f'\nOB02 user area {USER_AREA:.0f} m²: AD01 {AD01_COVER:.0%}, AD02 {AD02_COVER:.0%}, AD03 {AD03_COVER:.0%} below 31 °C')
    print('\nSensitivity (open sun / AD01):')
    for label, a, b in SENS:
        print(f'  {label:28s} {a:5.1f} / {b:5.1f} °C')
