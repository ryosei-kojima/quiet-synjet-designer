import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# --- ページ設定 ---
st.set_page_config(layout="wide", page_title="SynJet Designer v40")
st.title("🔇 Quiet SynJet Designer v40 (Trend Simulator)")

# --- 定数 ---
NU = 1.562e-5
C_SOUND = 343.0
C_CALIB = 23.9
POWER_U = 5.0
POWER_D = 10.0
BROADBAND_FLOOR = -15.0
BG_NOISE = 31.0
MEAS_DIST = 565.0
SP_NOISE_FACTOR = 38.0 # 固定

def get_a_weighting(f):
    if f <= 0: return -100.0
    f2 = f**2
    num = 12194**2 * f**4
    den = (f2 + 20.6**2) * np.sqrt((f2 + 107.7**2) * (f2 + 737.9**2)) * (f2 + 12194**2)
    gain = num / den
    val = 20 * np.log10(gain) + 2.0
    return val

# --- サイドバー ---
with st.sidebar:
    st.header("1. Settings (Fixed Factor=38)")
    
    st.header("2. Geometry Target")
    target_x = st.slider("Target Distance x (mm)", 50, 500, 150, 10)
    calc_D_mm = st.slider("Aperture Diameter D (mm)", 5.0, 200.0, 24.0, 1.0)
    
    xd_val = target_x / calc_D_mm
    st.metric("Ratio x/D", f"{xd_val:.1f}")
    if xd_val > 8.0: st.error("⚠️ x/D > 8")
    elif xd_val > 4.0: st.success("✅ Effective")
    else: st.info("⚡ Near Field")

    st.header("3. Speaker Spec")
    D_sp_mm = st.number_input("Effective Diameter (mm)", value=80.0)
    m0_g = st.number_input("Moving Mass m0 (g)", value=3.1)
    f0_hz = st.number_input("Resonance f0 (Hz)", value=82.0)
    Qts_fixed = st.number_input("Effective Qts", value=0.93, step=0.01)
    Z_ohm = st.number_input("Impedance (Ω)", value=8.0)
    
    st.header("4. Acoustics")
    h_cavity = st.number_input("Cavity Height (mm)", value=7.0)
    l_thickness = st.number_input("Plate Thickness (mm)", value=1.0)
    
    st.header("5. Drive")
    volts_pp = st.slider("Voltage (Vpp)", 1.0, 50.0, 5.0, 0.5)
    f_drive = st.slider("Frequency (Hz)", 10, 200, 63, 1)

# ==========================================================
# 物理計算
# ==========================================================
m0 = m0_g / 1000.0
Sd_m2 = (np.pi/4) * (D_sp_mm/1000.0)**2
Sd_mm2 = Sd_m2 * 1e6 
Sd_cm2 = Sd_m2 * 10000.0

w0 = 2 * np.pi * f0_hz
K_sp = m0 * w0**2 
Re_coil = Z_ohm * 0.8 
BL = np.sqrt( (w0 * m0 * Re_coil) / Qts_fixed )
F_peak = BL * ( (volts_pp/2) / Re_coil )

ratio = f_drive / f0_hz
mag_factor = 1 / np.sqrt( (1 - ratio**2)**2 + ( (1/Qts_fixed) * ratio )**2 )
X_static = F_peak / K_sp
X_peak_mm = (X_static * mag_factor) * 1000.0

Stroke_mm = 2.0 * X_peak_mm
V_slug_mm3 = Sd_mm2 * Stroke_mm 
V_slug_cm3 = V_slug_mm3 / 1000.0

L0_mm = (4 * V_slug_mm3) / (np.pi * calc_D_mm**2)
SR = L0_mm / calc_D_mm
U_exit = 2 * f_drive * (L0_mm / 1000.0) 
Re_val = (U_exit * (calc_D_mm/1000.0)) / NU

# ==========================================================
# 音響計算
# ==========================================================
# 1. Aero
St_aero = 0.2
f_noise = St_aero * U_exit / (calc_D_mm/1000.0)
raw_A_aero = get_a_weighting(f_noise)
eff_A_aero = max(raw_A_aero, BROADBAND_FLOOR)

if U_exit > 0.1:
    term_U = (10 * POWER_U) * np.log10(U_exit)
    term_D = POWER_D * np.log10(calc_D_mm/1000.0)
    term_R = 20 * np.log10(MEAS_DIST/1000.0)
    SPL_Z = C_CALIB + term_U + term_D - term_R
    L_aero_dBA = SPL_Z + eff_A_aero
else:
    L_aero_dBA = -99.9

# 2. Speaker
accel = (f_drive**2) * (X_peak_mm / 1000.0) 
L_sp_raw = 20 * np.log10(accel + 1e-9) + SP_NOISE_FACTOR 
A_drive = get_a_weighting(f_drive) 
L_sp_dBA = L_sp_raw + A_drive

# 3. Total
E_total = 10**(L_aero_dBA/10.0) + 10**(L_sp_dBA/10.0) + 10**(BG_NOISE/10.0)
Final_dBA = 10 * np.log10(E_total)

# Helmholtz
D_cav = D_sp_mm / 1000.0
Vc = (np.pi/4) * D_cav**2 * (h_cavity/1000.0)
A_port = (np.pi/4) * (calc_D_mm/1000.0)**2
Sigma = A_port / ((np.pi/4) * D_cav**2)
Le = (calc_D_mm/1000.0) * (Sigma**0.25) * (1.367 + 0.5173 * (l_thickness/calc_D_mm))
fH = (C_SOUND / (2*np.pi)) * np.sqrt(A_port / (Vc * Le))

detuning_ratio = f_drive / fH
if detuning_ratio < 0.95:
    res_gain = 1.0 / (1.0 - detuning_ratio**2)
    res_gain_db = 20 * np.log10(res_gain)
else:
    res_gain_db = 20.0

# ==========================================================
# 表示レイアウト
# ==========================================================

col1, col2, col3 = st.columns(3)

with col1:
    st.subheader("1. Design Check")
    if xd_val > 8.0: st.error(f"x/D = {xd_val:.1f}")
    else: st.success(f"x/D = {xd_val:.1f}")
        
    if SR < 1.0: st.error(f"SR = {SR:.2f}")
    else: st.metric("Stroke Ratio (SR)", f"{SR:.2f}", delta="OK")

with col2:
    st.subheader("2. Performance")
    st.metric("Exit Velocity", f"{U_exit:.2f} m/s")
    st.caption(f"Volume: {int(V_slug_mm3)} mm³")

with col3:
    st.subheader("3. Total Noise")
    if Final_dBA > 55: spl_col = "inverse"
    elif Final_dBA > 45: spl_col = "off"
    else: spl_col = "normal"
    st.metric("Est. SPL", f"{Final_dBA:.1f} dBA", delta_color=spl_col)
    
    st.markdown(f"""
    <small>
    Aero: {L_aero_dBA:.1f} dBA (f_n={int(f_noise)}Hz)<br>
    Spkr: {L_sp_dBA:.1f} dBA (f_d={int(f_drive)}Hz)<br>
    Back: {BG_NOISE} dBA
    </small>
    """, unsafe_allow_html=True)
    
    st.metric("Resonance Gain", f"+{res_gain_db:.1f} dB", delta=f"fH: {int(fH)}Hz")

st.divider()

# --- グラフ描画 ---
c1, c2 = st.columns(2)

with c1:
    fig, ax = plt.subplots(figsize=(6, 4))
    freqs = np.logspace(np.log10(10), np.log10(200), 200)
    rats = freqs / f0_hz
    mags = 1 / np.sqrt( (1 - rats**2)**2 + ( (1/Qts_fixed) * rats )**2 )
    
    Xs_static = (F_peak / K_sp) * 1000.0
    Vols_mm3 = Sd_mm2 * (2.0 * Xs_static * mags)
    
    ax.semilogx(freqs, Vols_mm3, 'b-', linewidth=2, label=f'Expelled Volume')
    ax.plot(f_drive, V_slug_mm3, 'ro', markersize=10, label='Current Point')
    
    req_vol = (np.pi/4 * calc_D_mm**3) 
    ax.axhline(y=req_vol, color='gray', linestyle='--', label='SR=1 Threshold')
    ax.axvline(x=f0_hz, color='k', linestyle=':', alpha=0.5)
    
    ax.set_xlabel('Frequency (Hz) [Log]')
    ax.set_ylabel('Expelled Volume (mm³)')
    ax.set_title(f'Speaker Response ({volts_pp}Vpp)')
    ax.set_xticks([10, 20, 50, 80, 100, 200])
    ax.set_xticklabels(['10', '20', '50', '80', '100', '200'])
    ax.grid(True, which="both", linestyle='--', alpha=0.5)
    ax.legend()
    st.pyplot(fig)

with c2:
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    
    # 背景ノイズレベルライン (基準としてこれだけ残す)
    ax2.axhline(y=BG_NOISE, color='gray', linestyle='--', linewidth=1.0, label='Background')

    # プロット (Total Noise)
    ax2.scatter(U_exit, Final_dBA, s=250, c='blue', edgecolors='black', zorder=10, label='Current Point')
    
    ax2.set_xlabel('Exit Velocity (m/s)')
    ax2.set_ylabel('Total Noise (dBA)')
    ax2.set_title('Velocity vs. Noise Level')
    
    ax2.set_xlim(0, 30)
    ax2.set_ylim(20, 60)
    
    ax2.grid(True, linestyle=':', alpha=0.7)
    ax2.legend(loc='lower right')
    
    st.pyplot(fig2)