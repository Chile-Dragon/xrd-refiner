"""
XRD 数据精修 — Web 全功能版
============================
功能对齐桌面版 xrd_refiner.py：
- 三种峰形: Gaussian / Lorentzian / PseudoVoigt
- 自动寻峰 / 手动添加删除 / 撤销重做
- 背景扣除 / 残差显示 / 单峰组分展示
- 实时预览 / 全谱拟合 / R² 评估
- 导出峰参数CSV / 精修数据TXT / 图谱PNG

运行: streamlit run xrd_web.py
"""

import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import signal, optimize
import copy, io, csv, re, base64
from io import BytesIO

# ============================================================
st.set_page_config(page_title="XRD 精修", page_icon="📊", layout="wide")

# ============================================================
#  核心类
# ============================================================
class DataManager:
    def __init__(self):
        self.x = self.y_raw = self.y_bg = self.y_corrected = None

    def load(self, file_bytes):
        text = file_bytes.decode('utf-8-sig', errors='replace')
        lines = [l.strip() for l in text.splitlines() if l.strip()]
        parsed = []
        for line in lines:
            nums = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', line)
            if len(nums) >= 2:
                parsed.append([float(nums[0]), float(nums[1])])
            elif len(nums) == 1:
                parsed.append([float(len(parsed)), float(nums[0])])
        if parsed:
            arr = np.array(parsed)
            self.x, self.y_raw = arr[:, 0], arr[:, 1]
        else:
            data = np.loadtxt(io.StringIO(text))
            self.x, self.y_raw = data[:, 0], data[:, 1]
        self.y_bg = np.zeros_like(self.y_raw)
        self.y_corrected = self.y_raw.copy()
        return len(self.x)

    def subtract_background(self, degree=5, window=80, n_iter=5):
        y_w = self.y_raw.copy()
        for _ in range(n_iter):
            nw = len(self.y_raw) // window
            xm, ym = [], []
            for i in range(nw):
                s, e = i * window, min((i + 1) * window, len(self.y_raw))
                idx = s + np.argmin(y_w[s:e])
                xm.append(self.x[idx]); ym.append(y_w[idx])
            xm, ym = np.array(xm), np.array(ym)
            bl = np.polyval(np.polyfit(xm, ym, degree), self.x)
            std = np.std(y_w - bl)
            y_w = np.where((y_w - bl) > 2.5 * std, bl + 2.5 * std, y_w)
        nw = len(self.y_raw) // window
        xm, ym = [], []
        for i in range(nw):
            s, e = i * window, min((i + 1) * window, len(self.y_raw))
            idx = s + np.argmin(y_w[s:e])
            xm.append(self.x[idx]); ym.append(y_w[idx])
        self.y_bg = np.polyval(np.polyfit(np.array(xm), np.array(ym), degree), self.x)
        self.y_corrected = np.maximum(self.y_raw - self.y_bg, 0)
        return self.y_bg, self.y_corrected

class PeakDetector:
    def __init__(self):
        self.peaks = []

    def auto_detect(self, x, y, prominence=None, distance=None):
        if prominence is None: prominence = np.max(y) * 0.03
        if distance is None: distance = max(3, len(x) // 300)
        idx_peaks, props = signal.find_peaks(y, prominence=prominence, distance=distance)
        heights = props.get('peak_heights', y[idx_peaks])
        self.peaks = []
        for k, i in enumerate(idx_peaks):
            h = float(heights[k]); half_h = h / 2
            left, right = i, i
            while left > 0 and y[left] > half_h: left -= 1
            while right < len(y) - 1 and y[right] > half_h: right += 1
            fwhm = x[right] - x[left]
            if fwhm <= 0: fwhm = 0.2
            self.peaks.append({'x': float(x[i]), 'height': h, 'fwhm': float(fwhm),
                               'area': h * fwhm * 1.064, 'type': 'Gaussian'})
        return self.peaks

    def add_peak(self, x_val, height=100):
        self.peaks.append({'x': float(x_val), 'height': float(height),
                           'fwhm': 0.3, 'area': height * 0.3 * 1.064, 'type': 'Gaussian'})
        self.peaks.sort(key=lambda p: p['x'])

    def remove_peak(self, index):
        if 0 <= index < len(self.peaks): del self.peaks[index]

    def find_nearest(self, x_val, x_tol=0.5):
        best_i, best_dist = -1, float('inf')
        for i, p in enumerate(self.peaks):
            d = abs(p['x'] - x_val)
            if d < best_dist: best_dist = d; best_i = i
        return (best_i, best_dist) if best_dist < x_tol else (-1, float('inf'))

    def update_peak(self, index, **kwargs):
        if 0 <= index < len(self.peaks):
            self.peaks[index].update(kwargs)
            self.peaks[index]['area'] = self.peaks[index]['height'] * self.peaks[index]['fwhm'] * 1.064

class PeakFitter:
    @staticmethod
    def gaussian(x, x0, A, sigma):
        return A * np.exp(-((x - x0) ** 2) / (2 * sigma ** 2))
    @staticmethod
    def lorentzian(x, x0, A, gamma):
        return A * gamma**2 / ((x - x0) ** 2 + gamma**2)
    @staticmethod
    def pseudo_voigt(x, x0, A, sigma, eta):
        g = A * np.exp(-((x - x0) ** 2) / (2 * sigma ** 2))
        l = A * sigma**2 / ((x - x0) ** 2 + sigma**2)
        return eta * l + (1 - eta) * g

    def fit(self, x, y, peaks, bg=None, method='Gaussian'):
        y_work = (y - bg) if bg is not None else y.copy()
        p0, lo, hi = [], [], []
        fn_map = {
            'Gaussian': (PeakFitter.gaussian, 3),
            'Lorentzian': (PeakFitter.lorentzian, 3),
            'PseudoVoigt': (PeakFitter.pseudo_voigt, 4),
        }
        fn_obj, n_p = fn_map[method]
        for p in peaks:
            s0 = p['fwhm'] / 2.355
            if method == 'PseudoVoigt':
                p0 += [p['x'], p['height'], max(s0, 0.02), 0.5]
                lo += [p['x'] - 1.0, 0, 0.005, 0.0]
                hi += [p['x'] + 1.0, p['height'] * 3, p['fwhm'] * 2, 1.0]
            else:
                p0 += [p['x'], p['height'], max(s0, 0.02)]
                lo += [p['x'] - 1.0, 0, 0.005]
                hi += [p['x'] + 1.0, p['height'] * 3, p['fwhm'] * 2]

        def multi(x_, *par):
            y_ = np.zeros_like(x_)
            for i in range(len(peaks)):
                y_ += fn_obj(x_, *par[i * n_p:(i + 1) * n_p])
            return y_

        try:
            popt, _ = optimize.curve_fit(multi, x, y_work, p0=p0, bounds=(lo, hi), maxfev=20000)
        except Exception:
            return None, None, None, 0

        y_fit_cmp = multi(x, *popt)
        y_fit_total = y_fit_cmp + (bg if bg is not None else 0)

        # 分解各峰组分
        components = []
        for i in range(len(peaks)):
            y_single = fn_obj(x, *popt[i * n_p:(i + 1) * n_p])
            if method == 'PseudoVoigt':
                x0, A, sigma, eta = popt[i * 4:(i + 1) * 4]
                fwhm = sigma * 2.355
            elif method == 'Lorentzian':
                x0, A, gamma = popt[i * 3:(i + 1) * 3]
                fwhm = gamma * 2.0
            else:  # Gaussian
                x0, A, sigma = popt[i * 3:(i + 1) * 3]
                fwhm = sigma * 2.355
            area = A * sigma * np.sqrt(2 * np.pi) if method != 'Lorentzian' else np.pi * A * gamma
            peaks[i].update({'x': float(x0), 'height': float(A),
                             'fwhm': float(fwhm), 'area': float(area), 'type': method})
            components.append(y_single)

        ss_res = np.sum((y_work - y_fit_cmp) ** 2)
        ss_tot = np.sum((y_work - np.mean(y_work)) ** 2)
        r = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        return y_fit_total, y_fit_cmp, components, r

def preview_fit(x, peaks, bg=None, method='Gaussian'):
    """快速预览（不优化）"""
    y = np.zeros_like(x)
    for p in peaks:
        s = p['fwhm'] / 2.355
        y += p['height'] * np.exp(-((x - p['x']) ** 2) / (2 * s ** 2))
    bg_arr = bg if bg is not None else np.zeros_like(x)
    return y + bg_arr

# ============================================================
#  Session State
# ============================================================
def init_session():
    defaults = {
        'data': None, 'detector': PeakDetector(), 'y_fit': None,
        'y_components': None, 'r_factor': 0.0, 'bg_enabled': False,
        'selected_peak': -1, 'undo_stack': [], 'redo_stack': [],
        'filename': None, 'show_residual': False, 'show_components': False,
        'fit_method': 'Gaussian',
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

def push_undo():
    st.session_state.undo_stack.append({
        'peaks': copy.deepcopy(st.session_state.detector.peaks),
        'sel': st.session_state.selected_peak,
        'bg': st.session_state.bg_enabled,
        'y_fit': copy.deepcopy(st.session_state.y_fit) if st.session_state.y_fit is not None else None,
        'y_components': copy.deepcopy(st.session_state.y_components) if st.session_state.y_components else None,
        'r': st.session_state.r_factor,
        'fit_method': st.session_state.fit_method,
    })
    if len(st.session_state.undo_stack) > 30:
        st.session_state.undo_stack.pop(0)
    st.session_state.redo_stack.clear()

def undo():
    if len(st.session_state.undo_stack) <= 1: return
    st.session_state.redo_stack.append(st.session_state.undo_stack.pop())
    restore_state(st.session_state.undo_stack[-1])

def redo():
    if not st.session_state.redo_stack: return
    s = st.session_state.redo_stack.pop()
    st.session_state.undo_stack.append(s)
    restore_state(s)

def restore_state(s):
    st.session_state.detector.peaks = copy.deepcopy(s['peaks'])
    st.session_state.selected_peak = s['sel']
    st.session_state.bg_enabled = s['bg']
    st.session_state.y_fit = copy.deepcopy(s['y_fit']) if s['y_fit'] else None
    st.session_state.y_components = copy.deepcopy(s['y_components']) if s['y_components'] else None
    st.session_state.r_factor = s['r']
    st.session_state.fit_method = s['fit_method']

# ============================================================
#  Plotly 图表
# ============================================================
COLORS_PEAK = [
    '#E6194B','#3CB44B','#FFE119','#4363D8','#F58231','#911EB4',
    '#46F0F0','#F032E6','#BCF60C','#FABEBE','#008080','#E6BEFF',
    '#9A6324','#FFFAC8','#800000','#AAFFC3','#808000','#FFD8B1',
    '#000075','#808080','#FFA500','#00CED1','#8B0000','#006400',
    '#FF69B4','#4B0082','#FFD700','#7CFC00','#00FFFF','#FF00FF',
]

def build_figure():
    dm = st.session_state.data
    if dm is None or dm.x is None:
        fig = go.Figure()
        fig.add_annotation(text="📂 请在侧边栏上传 XRD 数据文件", xref="paper", yref="paper",
                           x=0.5, y=0.5, showarrow=False, font=dict(size=20))
        fig.update_layout(height=650)
        return fig

    y_disp = dm.y_corrected if st.session_state.bg_enabled else dm.y_raw
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # 原始数据
    fig.add_trace(go.Scatter(x=dm.x, y=y_disp, mode='lines', name='XRD 数据',
        line=dict(color='#222222', width=0.8), hovertemplate='2θ=%{x:.3f}  I=%{y:.0f}<extra></extra>'),
        secondary_y=False)

    # 背景基线
    if st.session_state.bg_enabled:
        fig.add_trace(go.Scatter(x=dm.x, y=dm.y_bg, mode='lines', name='背景基线',
            line=dict(color='#CC4444', width=0.8, dash='dash')), secondary_y=False)

    # 拟合总曲线
    if st.session_state.y_fit is not None:
        fig.add_trace(go.Scatter(x=dm.x, y=st.session_state.y_fit, mode='lines', name='拟合曲线',
            line=dict(color='#0066CC', width=1.8)), secondary_y=False)

        # 各峰组分
        if st.session_state.show_components and st.session_state.y_components:
            for i, yc in enumerate(st.session_state.y_components):
                color = COLORS_PEAK[i % len(COLORS_PEAK)]
                fig.add_trace(go.Scatter(x=dm.x, y=yc, mode='lines',
                    name=f'峰{i+1}', line=dict(color=color, width=0.6, dash='dot'),
                    opacity=0.6), secondary_y=False)

    # 残差
    if st.session_state.show_residual and st.session_state.y_fit is not None:
        res = y_disp - st.session_state.y_fit
        fig.add_trace(go.Scatter(x=dm.x, y=res, mode='lines', name='残差',
            line=dict(color='#999999', width=0.5)), secondary_y=True)

    # 峰标记
    for i, p in enumerate(st.session_state.detector.peaks):
        sel = (i == st.session_state.selected_peak)
        color = '#E60000' if sel else '#E07000'
        size = 12 if sel else 7
        symbol = 'diamond' if sel else 'circle'
        fig.add_trace(go.Scatter(
            x=[p['x']], y=[p['height']], mode='markers+text',
            text=[f"{p['x']:.2f}"], textposition='top center',
            textfont=dict(size=8 if sel else 7, color=color, family='Arial'),
            marker=dict(color='#FF4400', size=size, symbol=symbol,
                        line=dict(color=color, width=2 if sel else 0.5)),
            name=f"峰{i+1}", showlegend=False,
            hovertemplate=f"<b>峰{i+1}</b><br>位置: {p['x']:.4f}°<br>强度: {p['height']:.1f}<br>FWHM: {p['fwhm']:.4f}<extra></extra>",
            customdata=[i],
        ), secondary_y=False)

        # 竖线
        fig.add_shape(type='line', x0=p['x'], x1=p['x'], y0=0, y1=p['height'],
                      line=dict(color=color, width=2 if sel else 0.6), layer='below')

    # 选中峰高度参考虚线
    if st.session_state.selected_peak >= 0 and st.session_state.selected_peak < len(st.session_state.detector.peaks):
        h = st.session_state.detector.peaks[st.session_state.selected_peak]['height']
        fig.add_hline(y=h, line=dict(color='#FF8800', width=1.2, dash='dot'),
                      annotation_text=f"高度={h:.0f}", annotation_position='top right')

    fig.update_layout(
        title=f"XRD 图谱 — {st.session_state.get('filename', '未加载')}",
        xaxis_title="2θ (度)", yaxis_title="强度 (counts)",
        height=650, hovermode='closest',
        legend=dict(x=0.99, y=0.99, xanchor='right', yanchor='top',
                    bgcolor='rgba(255,255,255,0.75)', font=dict(size=10)),
        margin=dict(l=50, r=30, t=50, b=50),
        clickmode='event+select',
    )
    fig.update_xaxes(rangeslider_visible=False)
    fig.update_yaxes(title_text="残差", secondary_y=True, showgrid=False)

    return fig

# ============================================================
#  UI
# ============================================================
st.title("📊 XRD 数据精修工具")
st.caption("自动寻峰 · 交互编辑 · 全谱拟合 · 数据导出  |  上传数据开始使用")
init_session()

# ── 侧边栏 ──
with st.sidebar:
    st.header("📂 数据加载")
    uploaded = st.file_uploader("拖放或选择 TXT/DAT/CSV 文件", type=['txt', 'dat', 'csv'],
                                 help="支持两列数据: 2θ角度 + 强度，Tab/逗号/空格分隔")
    if uploaded is not None:
        if (st.session_state.data is None or st.session_state.filename != uploaded.name):
            dm_new = DataManager()
            dm_new.load(uploaded.getvalue())
            st.session_state.data = dm_new
            st.session_state.detector = PeakDetector()
            st.session_state.y_fit = None; st.session_state.y_components = None; st.session_state.r_factor = 0
            st.session_state.selected_peak = -1
            st.session_state.filename = uploaded.name
            st.session_state.undo_stack = []; st.session_state.redo_stack = []
            push_undo()
            st.rerun()

    if st.session_state.data is not None and st.session_state.filename:
        dm = st.session_state.data
        st.success(f"✅ {st.session_state.filename}\n{dm.x[0]:.2f}° ~ {dm.x[-1]:.2f}° | {len(dm.x)} 点 | 强度 {dm.y_raw.min():.0f}~{dm.y_raw.max():.0f}")

    st.divider()

    # 寻峰
    st.header("🔍 寻峰参数")
    prominence_pct = st.slider("灵敏度 (%)", 0.5, 20.0, 3.0, 0.5, help="越小越灵敏, 找到更多峰")

    col_f1, col_f2 = st.columns(2)
    with col_f1:
        if st.button("🔍 自动寻峰", use_container_width=True, disabled=(st.session_state.data is None)):
            push_undo()
            dm = st.session_state.data
            y = dm.y_corrected if st.session_state.bg_enabled else dm.y_raw
            st.session_state.detector.auto_detect(
                dm.x, y, prominence=y.max() * prominence_pct / 100,
                distance=max(3, len(y) // 300))
            st.session_state.selected_peak = -1
            st.session_state.y_fit = None; st.session_state.y_components = None; st.session_state.r_factor = 0
            st.rerun()
    with col_f2:
        if st.button("📈 全谱拟合", use_container_width=True,
                      disabled=(st.session_state.data is None or not st.session_state.detector.peaks)):
            push_undo()
            dm = st.session_state.data
            y = dm.y_raw
            bg = dm.y_bg if st.session_state.bg_enabled else None
            yt, yc, comps, r = PeakFitter().fit(
                dm.x, y, st.session_state.detector.peaks, bg=bg,
                method=st.session_state.fit_method)
            if yt is not None:
                st.session_state.y_fit = yt
                st.session_state.y_components = comps
                st.session_state.r_factor = r
            else:
                st.error("拟合失败，请检查峰参数")
            st.rerun()

    st.divider()

    # 峰形选择
    st.header("⚙️ 拟合设置")
    method = st.selectbox("峰形函数", ['Gaussian', 'Lorentzian', 'PseudoVoigt'],
                           index=['Gaussian', 'Lorentzian', 'PseudoVoigt'].index(
                               st.session_state.fit_method))
    if method != st.session_state.fit_method:
        st.session_state.fit_method = method
        st.session_state.y_fit = None; st.session_state.y_components = None; st.session_state.r_factor = 0
        st.rerun()

    bg_on = st.checkbox("背景扣除", value=st.session_state.bg_enabled)
    if bg_on != st.session_state.bg_enabled:
        st.session_state.bg_enabled = bg_on
        if bg_on: st.session_state.data.subtract_background()
        st.session_state.y_fit = None; st.session_state.y_components = None; st.session_state.r_factor = 0
        st.rerun()

    resid_on = st.checkbox("显示残差", value=st.session_state.show_residual)
    if resid_on != st.session_state.show_residual:
        st.session_state.show_residual = resid_on; st.rerun()

    comps_on = st.checkbox("显示单峰组分", value=st.session_state.show_components)
    if comps_on != st.session_state.show_components:
        st.session_state.show_components = comps_on; st.rerun()

    st.divider()

    # 撤销/重做
    st.header("↩ 历史")
    col_u1, col_u2 = st.columns(2)
    with col_u1:
        st.button("↩ 撤销", use_container_width=True, on_click=undo,
                  disabled=len(st.session_state.undo_stack) <= 1)
    with col_u2:
        st.button("↪ 重做", use_container_width=True, on_click=redo,
                  disabled=not st.session_state.redo_stack)

    st.divider()

    # 导出
    st.header("📥 导出")
    if st.session_state.data is not None:
        dm = st.session_state.data

        # CSV
        if st.session_state.detector.peaks:
            csv_buf = io.StringIO()
            w = csv.writer(csv_buf)
            w.writerow(['#', '位置(°)', '强度', 'FWHM', '面积', '峰形'])
            for i, p in enumerate(st.session_state.detector.peaks):
                w.writerow([i+1, f"{p['x']:.4f}", f"{p['height']:.1f}",
                            f"{p['fwhm']:.4f}", f"{p['area']:.1f}", p.get('type','Gaussian')])
            st.download_button("📋 峰参数 (CSV)", csv_buf.getvalue(),
                               f"{st.session_state.filename or 'xrd'}_peaks.csv",
                               "text/csv", use_container_width=True)

        # TXT
        y_raw = dm.y_raw
        y_bg = dm.y_bg if st.session_state.bg_enabled else np.zeros_like(dm.x)
        y_corr = dm.y_corrected if st.session_state.bg_enabled else y_raw
        y_fit = st.session_state.y_fit if st.session_state.y_fit is not None else np.full_like(dm.x, np.nan)
        residual = y_corr - y_fit
        txt_buf = io.StringIO()
        txt_buf.write("2Theta(deg)\tRawIntensity\tBackground\tCorrected\tFitted\tResidual\n")
        for i in range(len(dm.x)):
            txt_buf.write(f"{dm.x[i]:.4f}\t{y_raw[i]:.2f}\t{y_bg[i]:.2f}\t{y_corr[i]:.2f}\t{y_fit[i]:.2f}\t{residual[i]:.2f}\n")
        st.download_button("📄 精修数据 (TXT)", txt_buf.getvalue(),
                           f"{st.session_state.filename or 'xrd'}_refined.txt",
                           "text/plain", use_container_width=True)

        # PNG
        try:
            fig_exp = build_figure()
            fig_exp.update_layout(width=1400, height=800)
            img_bytes = fig_exp.to_image(format='png', scale=2)
            st.download_button("🖼 图谱截图 (PNG)", img_bytes,
                               "xrd_chart.png", "image/png", use_container_width=True)
        except Exception:
            st.caption("PNG导出需要 kaleido: pip install kaleido")

# ── 主区域 ──
tab_chart, tab_table = st.tabs(["📊 图谱", "📋 峰列表"])

with tab_chart:
    fig = build_figure()
    chart_result = st.plotly_chart(fig, use_container_width=True, key="main_plot",
                                    on_select="rerun", selection_mode="points")

    # 处理图表点击
    if chart_result and chart_result.selection and chart_result.selection.points:
        pts = chart_result.selection.points
        if pts and pts[0].get('curve_number') is not None:
            # 找到被点击的峰标记
            pt = pts[0]
            click_x = pt.get('x', 0)
            # 在峰标记trace中查找（customdata包含峰索引）
            customdata = pt.get('customdata', None)
            if customdata is not None and isinstance(customdata, (int, float)):
                idx = int(customdata)
                if 0 <= idx < len(st.session_state.detector.peaks):
                    st.session_state.selected_peak = idx
                    st.rerun()

with tab_table:
    col_t1, col_t2, col_t3 = st.columns([1, 1, 1])
    with col_t1:
        st.metric("峰数", len(st.session_state.detector.peaks))
    with col_t2:
        st.metric("R²", f"{st.session_state.r_factor:.4f}" if st.session_state.r_factor > 0 else "-")
    with col_t3:
        st.metric("峰形", st.session_state.fit_method)

    # 峰列表
    if st.session_state.detector.peaks:
        peak_data = []
        for i, p in enumerate(st.session_state.detector.peaks):
            sel = "🔴" if i == st.session_state.selected_peak else ""
            peak_data.append({
                "": sel, "#": i+1, "位置(°)": f"{p['x']:.4f}",
                "强度": f"{p['height']:.1f}", "FWHM": f"{p['fwhm']:.4f}",
                "面积": f"{p['area']:.1f}",
            })
        st.dataframe(peak_data, use_container_width=True, hide_index=True,
                      column_config={"": st.column_config.TextColumn("选", width="small")})
    else:
        st.info("暂无峰数据，请先自动寻峰或手动添加。")

# ── 编辑面板 ──
st.divider()
st.subheader("✏️ 编辑峰" if st.session_state.selected_peak >= 0 else "✏️ 操作")

if st.session_state.selected_peak >= 0 and st.session_state.selected_peak < len(st.session_state.detector.peaks):
    p = st.session_state.detector.peaks[st.session_state.selected_peak]
    dm = st.session_state.data
    idx = st.session_state.selected_peak

    ed1, ed2, ed3, ed4 = st.columns([2, 2, 2, 1])

    with ed1:
        new_x = st.number_input("位置 (°)", value=float(p['x']),
            step=0.01, min_value=float(dm.x[0]), max_value=float(dm.x[-1]),
            format="%.4f", key="edit_x")

    with ed2:
        new_h = st.number_input("强度", value=float(p['height']),
            step=10.0, min_value=0.0, format="%.1f", key="edit_h")

    with ed3:
        new_w = st.number_input("FWHM", value=float(p['fwhm']),
            step=0.01, min_value=0.005, format="%.4f", key="edit_w")

    with ed4:
        st.write("")  # spacer
        st.write("")
        if st.button("✅ 应用", use_container_width=True, type="primary"):
            push_undo()
            st.session_state.detector.update_peak(
                idx, x=new_x, height=new_h, fwhm=new_w)
            dm = st.session_state.data
            st.session_state.y_fit = preview_fit(
                dm.x, st.session_state.detector.peaks,
                bg=dm.y_bg if st.session_state.bg_enabled else None)
            st.session_state.y_components = None; st.session_state.r_factor = 0
            st.rerun()

    # 快捷微调按钮
    btn1, btn2, btn3, btn4, btn5 = st.columns(5)
    with btn1:
        if st.button("◀◀ -1°", use_container_width=True, key="nx_big"):
            push_undo()
            st.session_state.detector.update_peak(idx, x=max(dm.x[0], p['x'] - 1))
            st.rerun()
    with btn2:
        if st.button("◀ -0.1°", use_container_width=True, key="nx_mid"):
            push_undo()
            st.session_state.detector.update_peak(idx, x=max(dm.x[0], p['x'] - 0.1))
            st.rerun()
    with btn3:
        if st.button("▶ +0.1°", use_container_width=True, key="nx_mid2"):
            push_undo()
            st.session_state.detector.update_peak(idx, x=min(dm.x[-1], p['x'] + 0.1))
            st.rerun()
    with btn4:
        if st.button("▶▶ +1°", use_container_width=True, key="nx_big2"):
            push_undo()
            st.session_state.detector.update_peak(idx, x=min(dm.x[-1], p['x'] + 1))
            st.rerun()
    with btn5:
        if st.button("🗑 删除", use_container_width=True, key="del_peak"):
            push_undo()
            st.session_state.detector.remove_peak(idx)
            st.session_state.selected_peak = -1
            st.session_state.y_fit = None; st.session_state.y_components = None
            st.rerun()

    # 强度微调
    bh1, bh2, bh3 = st.columns(3)
    with bh1:
        if st.button("🔽 强度 -100", use_container_width=True):
            push_undo()
            st.session_state.detector.update_peak(idx, height=max(1, p['height'] - 100))
            st.rerun()
    with bh2:
        if st.button("🔼 强度 +100", use_container_width=True):
            push_undo()
            st.session_state.detector.update_peak(idx, height=p['height'] + 100)
            st.rerun()
    with bh3:
        fwhm_adj = st.number_input("FWHM缩放", value=1.0, step=0.1, min_value=0.1, max_value=5.0, key="fwhm_scale")
        if st.button("应用FWHM", use_container_width=True):
            push_undo()
            st.session_state.detector.update_peak(idx, fwhm=p['fwhm'] * fwhm_adj)
            st.rerun()

else:
    st.info("👆 点击图谱上的峰标记 / 右侧峰列表来选择峰，然后在此处编辑参数。")

# 添加峰
st.divider()
st.subheader("➕ 添加/删除峰")
if st.session_state.data is not None:
    dm = st.session_state.data
    col_a1, col_a2, col_a3 = st.columns([2, 1, 1])
    with col_a1:
        add_x = st.number_input("新峰位置 (°)", value=float((dm.x[0]+dm.x[-1])/2),
            min_value=float(dm.x[0]), max_value=float(dm.x[-1]), format="%.4f", key="add_x")
    with col_a2:
        if st.button("➕ 添加峰", use_container_width=True, type="primary"):
            push_undo()
            idx = np.argmin(np.abs(dm.x - add_x))
            y_val = dm.y_corrected[idx] if st.session_state.bg_enabled else dm.y_raw[idx]
            st.session_state.detector.add_peak(add_x, max(y_val, 1))
            st.session_state.selected_peak = len(st.session_state.detector.peaks) - 1
            st.session_state.y_fit = preview_fit(
                dm.x, st.session_state.detector.peaks,
                bg=dm.y_bg if st.session_state.bg_enabled else None)
            st.session_state.y_components = None; st.session_state.r_factor = 0
            st.rerun()
    with col_a3:
        if st.button("🗑 清空全部峰", use_container_width=True):
            push_undo()
            st.session_state.detector.peaks.clear()
            st.session_state.selected_peak = -1
            st.session_state.y_fit = None; st.session_state.y_components = None; st.session_state.r_factor = 0
            st.rerun()

# 底部状态栏
if st.session_state.data is not None:
    st.divider()
    n = len(st.session_state.data.x)
    npk = len(st.session_state.detector.peaks)
    rs = f" | R²={st.session_state.r_factor:.4f}" if st.session_state.r_factor > 0 else ""
    bgs = " | 背景ON" if st.session_state.bg_enabled else ""
    comps_s = " | 显示组分" if st.session_state.show_components else ""
    st.caption(f"📊 {st.session_state.filename}: {n} 点, 2θ={st.session_state.data.x[0]:.2f}°~{st.session_state.data.x[-1]:.2f}° "
               f"| 🔴 {npk} 峰{rs}{bgs}{comps_s} "
               f"| 峰形: {st.session_state.fit_method}")
