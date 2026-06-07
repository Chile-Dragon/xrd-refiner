"""
XRD 数据精修 — Web 版
=====================
基于 Streamlit + Plotly，在线进行 XRD 图谱寻峰、拟合、精修。
运行: streamlit run xrd_web.py
"""

import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import signal, optimize
import copy, io, csv

# ============================================================
#  核心类 (从桌面版复用)
# ============================================================
class DataManager:
    def __init__(self):
        self.x = self.y_raw = self.y_bg = self.y_corrected = None

    def load(self, file_bytes, filename):
        import re
        text = file_bytes.decode('utf-8-sig', errors='replace')
        lines = [l.strip() for l in text.splitlines() if l.strip()]
        parsed = []
        for line in lines:
            nums = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', line)
            if len(nums) >= 2:
                parsed.append([float(nums[0]), float(nums[1])])
            elif len(nums) == 1:
                parsed.append([float(len(parsed)), float(nums[0])])
        if not parsed:
            data = np.loadtxt(io.StringIO(text))
            self.x, self.y_raw = data[:, 0], data[:, 1]
        else:
            arr = np.array(parsed)
            self.x, self.y_raw = arr[:, 0], arr[:, 1]
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

    def fit(self, x, y, peaks, bg=None, method='Gaussian'):
        y_work = (y - bg) if bg is not None else y.copy()
        p0, lo, hi = [], [], []
        n_p = 3  # Gaussian only for simplicity
        for p in peaks:
            s0 = p['fwhm'] / 2.355
            p0 += [p['x'], p['height'], max(s0, 0.02)]
            lo += [p['x'] - 1.0, 0, 0.005]
            hi += [p['x'] + 1.0, p['height'] * 3, p['fwhm'] * 2]

        def multi(x_, *par):
            y_ = np.zeros_like(x_)
            for i in range(len(peaks)):
                y_ += PeakFitter.gaussian(x_, *par[i * 3:(i + 1) * 3])
            return y_

        try:
            popt, _ = optimize.curve_fit(multi, x, y_work, p0=p0, bounds=(lo, hi), maxfev=20000)
        except Exception:
            return None, None, 0

        y_fit_cmp = multi(x, *popt)
        y_fit_total = y_fit_cmp + (bg if bg is not None else 0)
        for i in range(len(peaks)):
            x0, A, sigma = popt[i * 3:(i + 1) * 3]
            fwhm = sigma * 2.355
            area = A * sigma * np.sqrt(2 * np.pi)
            peaks[i].update({'x': float(x0), 'height': float(A),
                             'fwhm': float(fwhm), 'area': float(area)})
        ss_res = np.sum((y_work - y_fit_cmp) ** 2)
        ss_tot = np.sum((y_work - np.mean(y_work)) ** 2)
        r = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        return y_fit_total, y_fit_cmp, r

def preview_fit(x, peaks, bg=None):
    """快速生成预览曲线"""
    y = np.zeros_like(x)
    for p in peaks:
        sigma = p['fwhm'] / 2.355
        y += p['height'] * np.exp(-((x - p['x']) ** 2) / (2 * sigma ** 2))
    return y + (bg if bg is not None else 0)

# ============================================================
#  Session State 初始化
# ============================================================
def init_session():
    defaults = {
        'data': None, 'detector': PeakDetector(), 'y_fit': None,
        'y_fit_components': None, 'r_factor': 0.0, 'bg_enabled': False,
        'selected_peak': -1, 'undo_stack': [], 'redo_stack': [],
        'filename': None, 'show_residual': False,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

def push_undo():
    st.session_state.undo_stack.append({
        'peaks': copy.deepcopy(st.session_state.detector.peaks),
        'sel': st.session_state.selected_peak,
        'bg': st.session_state.bg_enabled,
        'y_fit': st.session_state.y_fit.copy() if st.session_state.y_fit is not None else None,
        'r': st.session_state.r_factor,
    })
    if len(st.session_state.undo_stack) > 30:
        st.session_state.undo_stack.pop(0)
    st.session_state.redo_stack.clear()

def undo():
    if len(st.session_state.undo_stack) <= 1: return
    st.session_state.redo_stack.append(st.session_state.undo_stack.pop())
    s = st.session_state.undo_stack[-1]
    restore_state(s)

def redo():
    if not st.session_state.redo_stack: return
    s = st.session_state.redo_stack.pop()
    st.session_state.undo_stack.append(s)
    restore_state(s)

def restore_state(s):
    st.session_state.detector.peaks = copy.deepcopy(s['peaks'])
    st.session_state.selected_peak = s['sel']
    st.session_state.bg_enabled = s['bg']
    st.session_state.y_fit = s['y_fit'].copy() if s['y_fit'] else None
    st.session_state.r_factor = s['r']

# ============================================================
#  Plotly 绘图
# ============================================================
def build_figure():
    dm = st.session_state.data
    if dm is None or dm.x is None:
        fig = go.Figure()
        fig.add_annotation(text="请上传 XRD 数据文件", xref="paper", yref="paper",
                           x=0.5, y=0.5, showarrow=False, font=dict(size=20))
        fig.update_layout(height=600)
        return fig

    y_disp = dm.y_corrected if st.session_state.bg_enabled else dm.y_raw
    fig = go.Figure()

    # 原始数据
    fig.add_trace(go.Scatter(
        x=dm.x, y=y_disp, mode='lines', name='XRD 数据',
        line=dict(color='#222222', width=0.8), hoverinfo='x+y',
    ))

    # 背景基线
    if st.session_state.bg_enabled:
        fig.add_trace(go.Scatter(
            x=dm.x, y=dm.y_bg, mode='lines', name='背景基线',
            line=dict(color='#CC4444', width=0.8, dash='dash'),
        ))

    # 拟合曲线
    if st.session_state.y_fit is not None:
        fig.add_trace(go.Scatter(
            x=dm.x, y=st.session_state.y_fit, mode='lines', name='拟合曲线',
            line=dict(color='#0066CC', width=1.5),
        ))

    # 残差
    if st.session_state.show_residual and st.session_state.y_fit is not None:
        res = y_disp - st.session_state.y_fit
        fig.add_trace(go.Scatter(
            x=dm.x, y=res + 5, mode='lines', name='残差',
            line=dict(color='#999999', width=0.5), yaxis='y2',
        ))

    # 峰标记
    for i, p in enumerate(st.session_state.detector.peaks):
        sel = (i == st.session_state.selected_peak)
        color = '#E60000' if sel else '#E07000'
        size = 10 if sel else 6
        symbol = 'diamond' if sel else 'circle'
        h_actual = p['height']
        fig.add_trace(go.Scatter(
            x=[p['x']], y=[h_actual], mode='markers+text',
            text=[f"{p['x']:.2f}"], textposition='top center',
            textfont=dict(size=7, color=color),
            marker=dict(color='#FF4400', size=size, symbol=symbol,
                        line=dict(color=color, width=1.5 if sel else 0.5)),
            name=f"峰{i+1}", showlegend=False,
            hovertemplate=f"峰{i+1}<br>位置: {p['x']:.3f}°<br>强度: {p['height']:.0f}<br>FWHM: {p['fwhm']:.3f}<extra></extra>",
        ))
        # 竖线
        fig.add_shape(type='line', x0=p['x'], x1=p['x'], y0=0, y1=p['height'],
                      line=dict(color=color, width=1.5 if sel else 0.6), layer='below')

    # 选中峰高度参考线
    if st.session_state.selected_peak >= 0 and st.session_state.selected_peak < len(st.session_state.detector.peaks):
        h = st.session_state.detector.peaks[st.session_state.selected_peak]['height']
        fig.add_hline(y=h, line=dict(color='#FF8800', width=1, dash='dot'), layer='below')

    fig.update_layout(
        title=f"XRD 图谱 — {st.session_state.get('filename', '未加载')}",
        xaxis_title="2θ (度)", yaxis_title="强度 (counts)",
        height=600, hovermode='closest',
        legend=dict(x=0.98, y=0.98, xanchor='right', yanchor='top',
                    bgcolor='rgba(255,255,255,0.7)'),
        margin=dict(l=50, r=20, t=50, b=50),
    )
    fig.update_xaxes(rangeslider_visible=False)

    if st.session_state.show_residual and st.session_state.y_fit is not None:
        fig.update_layout(yaxis2=dict(overlaying='y', side='right', range=[-10, 50]))

    return fig

# ============================================================
#  UI
# ============================================================
st.set_page_config(page_title="XRD 数据精修", page_icon="📊", layout="wide")
st.title("📊 XRD 数据精修工具 (Web)")
init_session()

# ── 侧边栏 ──
with st.sidebar:
    st.header("📂 数据")
    uploaded = st.file_uploader("上传 XRD 数据文件 (.txt)", type=['txt', 'dat', 'csv'])
    if uploaded is not None:
        if st.session_state.data is None or st.session_state.filename != uploaded.name:
            dm_new = DataManager()
            dm_new.load(uploaded.getvalue(), uploaded.name)
            st.session_state.data = dm_new
            st.session_state.detector = PeakDetector()
            st.session_state.y_fit = None; st.session_state.r_factor = 0
            st.session_state.selected_peak = -1
            st.session_state.filename = uploaded.name
            st.session_state.undo_stack = []; st.session_state.redo_stack = []
            push_undo()
            st.success(f"已加载: {uploaded.name}")
            st.rerun()

    if st.session_state.data is not None:
        st.info(f"📄 {st.session_state.filename}\n{len(st.session_state.data.x)} 个数据点")

    st.header("🔍 寻峰")
    prominence_pct = st.slider("灵敏度 (%)", 0.5, 20.0, 3.0, 0.5)
    col1, col2 = st.columns(2)
    with col1:
        if st.button("🔍 自动寻峰", use_container_width=True):
            push_undo()
            y = st.session_state.data.y_corrected if st.session_state.bg_enabled else st.session_state.data.y_raw
            st.session_state.detector.auto_detect(
                st.session_state.data.x, y,
                prominence=y.max() * prominence_pct / 100,
                distance=max(3, len(y) // 300))
            st.session_state.selected_peak = -1
            st.session_state.y_fit = None; st.session_state.r_factor = 0
            st.rerun()
    with col2:
        if st.button("📈 全谱拟合", use_container_width=True):
            push_undo()
            y = st.session_state.data.y_raw
            bg = st.session_state.data.y_bg if st.session_state.bg_enabled else None
            yt, yc, r = PeakFitter().fit(
                st.session_state.data.x, y, st.session_state.detector.peaks,
                bg=bg, method='Gaussian')
            if yt is not None:
                st.session_state.y_fit = yt
                st.session_state.y_fit_components = yc
                st.session_state.r_factor = r
                st.success(f"拟合完成: R²={r:.4f}")
            else:
                st.error("拟合失败")
            st.rerun()

    st.header("⚙️ 设置")
    bg_on = st.checkbox("背景扣除", value=st.session_state.bg_enabled)
    if bg_on != st.session_state.bg_enabled:
        st.session_state.bg_enabled = bg_on
        if bg_on: st.session_state.data.subtract_background()
        push_undo()
        st.session_state.y_fit = None; st.session_state.r_factor = 0
        st.rerun()

    resid_on = st.checkbox("显示残差", value=st.session_state.show_residual)
    if resid_on != st.session_state.show_residual:
        st.session_state.show_residual = resid_on
        st.rerun()

    st.header("↩ 操作")
    c1, c2 = st.columns(2)
    with c1:
        if st.button("↩ 撤销", use_container_width=True): undo(); st.rerun()
    with c2:
        if st.button("↪ 重做", use_container_width=True): redo(); st.rerun()

    if st.button("🗑 删除选中峰", use_container_width=True):
        if st.session_state.selected_peak >= 0:
            push_undo()
            st.session_state.detector.remove_peak(st.session_state.selected_peak)
            st.session_state.selected_peak = -1
            st.session_state.y_fit = None
            st.rerun()

    st.header("📥 导出")
    if st.session_state.data is not None:
        # CSV 峰参数
        if st.session_state.detector.peaks:
            csv_buf = io.StringIO()
            w = csv.writer(csv_buf)
            w.writerow(['#', '位置(°)', '强度', 'FWHM', '面积', '峰形'])
            for i, p in enumerate(st.session_state.detector.peaks):
                w.writerow([i+1, f"{p['x']:.4f}", f"{p['height']:.1f}",
                            f"{p['fwhm']:.4f}", f"{p['area']:.1f}", p['type']])
            st.download_button("📋 峰参数 (CSV)", csv_buf.getvalue(),
                               "xrd_peaks.csv", "text/csv", use_container_width=True)

        # TXT 精修数据
        dm = st.session_state.data
        y_raw = dm.y_raw
        y_bg = dm.y_bg if st.session_state.bg_enabled else np.zeros_like(dm.x)
        y_corr = dm.y_corrected if st.session_state.bg_enabled else y_raw
        y_fit = st.session_state.y_fit if st.session_state.y_fit is not None else np.full_like(dm.x, np.nan)
        residual = y_corr - y_fit

        txt_buf = io.StringIO()
        txt_buf.write("2Theta(deg)\tRawIntensity\tBackground\tCorrected\tFitted\tResidual\n")
        for i in range(len(dm.x)):
            txt_buf.write(f"{dm.x[i]:.4f}\t{y_raw[i]:.4f}\t{y_bg[i]:.4f}\t{y_corr[i]:.4f}\t{y_fit[i]:.4f}\t{residual[i]:.4f}\n")
        st.download_button("📄 精修数据 (TXT)", txt_buf.getvalue(),
                           "xrd_refined_data.txt", "text/plain", use_container_width=True)

# ── 主区域 ──
col_chart, col_edit = st.columns([3, 1])

with col_chart:
    fig = build_figure()

    # 点击事件：选中峰
    clicked = st.plotly_chart(fig, use_container_width=True, key="main_chart",
                               on_select="rerun", selection_mode="points")

    if clicked and clicked.selection:
        pts = clicked.selection.points
        if pts:
            # 找到最近的峰
            click_x = pts[0].get('x', 0)
            idx, _ = st.session_state.detector.find_nearest(click_x, x_tol=1.0)
            if idx >= 0:
                st.session_state.selected_peak = idx
                st.rerun()

with col_edit:
    st.subheader("📝 峰列表")

    # 峰数量统计
    n_peaks = len(st.session_state.detector.peaks)
    r_str = f"R²={st.session_state.r_factor:.4f}" if st.session_state.r_factor > 0 else ""
    st.caption(f"{n_peaks} 个峰 | {r_str}")

    # 快捷跳转按钮
    for i, p in enumerate(st.session_state.detector.peaks):
        sel = (i == st.session_state.selected_peak)
        btn_label = f"{'🔴' if sel else '  '} #{i+1}: {p['x']:.2f}° | {p['height']:.0f}"
        if st.button(btn_label, use_container_width=True, key=f"sel_{i}",
                     type="primary" if sel else "secondary"):
            st.session_state.selected_peak = i
            st.rerun()

    st.divider()

    # 选中峰编辑
    if st.session_state.selected_peak >= 0 and st.session_state.selected_peak < n_peaks:
        st.subheader("✏️ 编辑峰")
        p = st.session_state.detector.peaks[st.session_state.selected_peak]
        dm = st.session_state.data

        new_x = st.number_input("位置 (°)", value=float(p['x']), step=0.01,
                                 min_value=float(dm.x[0]), max_value=float(dm.x[-1]),
                                 format="%.3f", key=f"edit_x")
        new_h = st.number_input("强度", value=float(p['height']), step=1.0,
                                 min_value=1.0, format="%.0f", key=f"edit_h")
        new_w = st.number_input("FWHM", value=float(p['fwhm']), step=0.01,
                                 min_value=0.01, format="%.3f", key=f"edit_w")

        col_a, col_b = st.columns(2)
        with col_a:
            if st.button("✅ 应用", use_container_width=True):
                push_undo()
                st.session_state.detector.update_peak(
                    st.session_state.selected_peak, x=new_x, height=new_h, fwhm=new_w)
                st.session_state.y_fit = preview_fit(
                    st.session_state.data.x, st.session_state.detector.peaks,
                    bg=st.session_state.data.y_bg if st.session_state.bg_enabled else None)
                st.rerun()
        with col_b:
            if st.button("🗑 删除", use_container_width=True):
                push_undo()
                st.session_state.detector.remove_peak(st.session_state.selected_peak)
                st.session_state.selected_peak = -1
                st.session_state.y_fit = None
                st.rerun()

        # 添加峰
        st.divider()
        st.subheader("➕ 添加峰")
        add_x = st.number_input("新峰位置 (°)", value=float(dm.x[len(dm.x)//2]),
                                 min_value=float(dm.x[0]), max_value=float(dm.x[-1]), format="%.3f")
        if st.button("添加", use_container_width=True):
            push_undo()
            idx = np.argmin(np.abs(dm.x - add_x))
            y_val = dm.y_corrected[idx] if st.session_state.bg_enabled else dm.y_raw[idx]
            st.session_state.detector.add_peak(add_x, max(y_val, 1))
            st.session_state.y_fit = preview_fit(
                dm.x, st.session_state.detector.peaks,
                bg=dm.y_bg if st.session_state.bg_enabled else None)
            st.rerun()

# ── 底部状态 ──
if st.session_state.data is not None:
    st.divider()
    n = len(st.session_state.data.x)
    npk = len(st.session_state.detector.peaks)
    rs = f" | R²={st.session_state.r_factor:.4f}" if st.session_state.r_factor > 0 else ""
    bgs = " | 背景ON" if st.session_state.bg_enabled else ""
    st.caption(f"数据: {n} 点 | 峰: {npk} 个{rs}{bgs} | 点击图谱上的峰标记即可选中编辑")
