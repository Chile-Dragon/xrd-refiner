#!/usr/bin/env python3
"""
XRD 数据精修可视化工具
======================
- 拖拽峰位 (水平移动)  |  Shift+拖拽 (调整强度)
- 右键添加峰  |  右键峰上→删除  |  双击表格直接编辑参数
- 多峰拟合 (Gaussian / Lorentzian / PseudoVoigt)
- 背景扣除、撤销/重做、CSV/PNG 导出
- 支持拖放 TXT 文件到窗口直接打开

依赖: numpy, scipy, matplotlib, tkinter
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
from scipy import signal, optimize
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib
import copy, os, csv, sys

matplotlib.use('TkAgg')

# ── 中文字体 ──
for fname in ['Microsoft YaHei', 'SimHei', 'SimSun', 'KaiTi', 'FangSong']:
    try:
        matplotlib.rcParams['font.sans-serif'] = [fname] + matplotlib.rcParams['font.sans-serif']
        matplotlib.rcParams['axes.unicode_minus'] = False
        break
    except:
        pass

# ============================================================
COLORS = {
    'raw':        '#222222',
    'bg':         '#CC4444',
    'fit':        '#0066CC',
    'residual':   '#999999',
    'peak_line':  '#E07000',
    'peak_dot':   '#FF4400',
    'selected':   '#E60000',
    'height_line':'#FF8800',
}
# ============================================================

class DataManager:
    def __init__(self):
        self.x = self.y_raw = self.y_bg = self.y_corrected = None
        self.filepath = None

    def load(self, filepath):
        """加载数据 — 自动处理多种格式（跳过非数字行）"""
        import warnings, re

        # 先读取原始内容判断格式
        with open(filepath, 'r', encoding='utf-8-sig', errors='replace') as f:
            lines = [l.strip() for l in f if l.strip()]

        # 手动解析：逐行提取数字对
        parsed = []
        for line in lines:
            # 跳过纯注释/标题行（不包含任何数字）
            if not re.search(r'[\d.]+', line):
                continue
            # 提取所有数字（支持科学计数法、负数）
            nums = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', line)
            if len(nums) >= 2:
                # 有至少两列 → 取前两列
                parsed.append([float(nums[0]), float(nums[1])])
            elif len(nums) == 1:
                # 单列 → 自动生成 x
                parsed.append([float(len(parsed)), float(nums[0])])

        if parsed:
            data = np.array(parsed)
        else:
            # 回退到 np.loadtxt
            data = np.loadtxt(filepath)

        self.x = data[:, 0]
        self.y_raw = data[:, 1]
        self.filepath = filepath
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
        self.peaks = []  # [{x, height, fwhm, area, type}, ...]

    def auto_detect(self, x, y, prominence=None, distance=None):
        if prominence is None: prominence = np.max(y) * 0.03
        if distance is None:   distance = max(3, len(x) // 300)
        idx_peaks, props = signal.find_peaks(y, prominence=prominence, distance=distance)
        heights = props.get('peak_heights', y[idx_peaks])
        self.peaks = []
        for k, i in enumerate(idx_peaks):
            h = float(heights[k])
            half_h = h / 2
            left = i
            while left > 0 and y[left] > half_h: left -= 1
            right = i
            while right < len(y) - 1 and y[right] > half_h: right += 1
            fwhm = x[right] - x[left]
            if fwhm <= 0: fwhm = 0.2
            self.peaks.append({
                'x': float(x[i]), 'height': h,
                'fwhm': float(fwhm), 'area': h * fwhm * 1.064,
                'type': 'Gaussian'
            })
        return self.peaks

    def add_peak(self, x_val, height=100):
        self.peaks.append({
            'x': float(x_val), 'height': float(height),
            'fwhm': 0.3, 'area': float(height) * 0.3 * 1.064, 'type': 'Gaussian'
        })
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
        fn = {'Gaussian': (PeakFitter.gaussian, 3),
              'Lorentzian': (PeakFitter.lorentzian, 3),
              'PseudoVoigt': (PeakFitter.pseudo_voigt, 4)}[method]
        fn_obj, n_p = fn
        for p in peaks:
            s0 = p['fwhm'] / 2.355
            if method == 'PseudoVoigt':
                p0 += [p['x'], p['height'], max(s0, 0.02), 0.5]
                lo += [p['x']-1.0, 0, 0.005, 0.0]
                hi += [p['x']+1.0, p['height']*3, p['fwhm']*2, 1.0]
            else:
                p0 += [p['x'], p['height'], max(s0, 0.02)]
                lo += [p['x']-1.0, 0, 0.005]
                hi += [p['x']+1.0, p['height']*3, p['fwhm']*2]
        def multi(x_, *par):
            y_ = np.zeros_like(x_)
            for i in range(len(peaks)): y_ += fn_obj(x_, *par[i*n_p:(i+1)*n_p])
            return y_
        try:
            popt, _ = optimize.curve_fit(multi, x, y_work, p0=p0, bounds=(lo, hi), maxfev=20000)
        except Exception as e:
            print(f"Fitting failed: {e}")
            return None, None, None, 0
        y_fit_cmp = multi(x, *popt)
        y_fit_total = y_fit_cmp + bg if bg is not None else y_fit_cmp
        for i in range(len(peaks)):
            if method == 'PseudoVoigt':
                x0, A, sigma, eta = popt[i*4:(i+1)*4]
                fwhm = sigma * 2.355
                area = A * sigma * np.sqrt(2*np.pi)
            else:
                x0, A, s = popt[i*3:(i+1)*3]
                fwhm = s*2.0 if method == 'Lorentzian' else s*2.355
                area = np.pi*A*s if method == 'Lorentzian' else A*s*np.sqrt(2*np.pi)
            peaks[i].update({'x': float(x0), 'height': float(A), 'fwhm': float(fwhm), 'area': float(area)})
        ss_res = np.sum((y_work - y_fit_cmp)**2)
        ss_tot = np.sum((y_work - np.mean(y_work))**2)
        r = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        return y_fit_total, y_fit_cmp, None, r

# ============================================================
class XrdRefinerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("XRD 数据精修工具")
        self.root.geometry("1400x850")

        self.data = DataManager()
        self.detector = PeakDetector()
        self.fitter = PeakFitter()

        self.selected_peak = -1
        self.dragging_peak = -1
        self.drag_mode = None      # 'x' or 'height'
        self.drag_start = None
        self.y_fit = None
        self.r_factor = 0.0
        self.bg_enabled = False
        self.show_residual = False
        self.fit_method = tk.StringVar(value='Gaussian')

        self.undo_stack = []; self.redo_stack = []; self._push_undo()

        self.raw_line = self.bg_line = self.fit_line = self.residual_line = None
        self.peak_lines = []; self.peak_dots = []; self.peak_labels = []
        self.height_line = None
        self.edit_entry = None  # 表格内联编辑控件

        self._build_ui()
        self._connect_events()
        self._setup_drag_drop()

    # ── UI ────────────────────────────────
    def _build_ui(self):
        main_pw = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_pw.pack(fill=tk.BOTH, expand=True)

        # 左侧
        left_frame = ttk.Frame(main_pw); main_pw.add(left_frame, weight=3)
        self.fig = Figure(figsize=(10, 7), dpi=100)
        self.ax_main = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=left_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        tb_frame = ttk.Frame(left_frame); tb_frame.pack(fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas, tb_frame); self.toolbar.update()

        # 右侧
        right_frame = ttk.Frame(main_pw, width=370); main_pw.add(right_frame, weight=1)

        # 文件信息
        self.file_label = ttk.Label(right_frame, text="拖放 TXT 文件到窗口 / 点击下方按钮加载",
                                     foreground='gray', font=('', 9))
        self.file_label.pack(pady=(8,2), padx=5, anchor='w')

        # 操作指南
        self.help_frame = ttk.LabelFrame(right_frame, text="📖 操作指南", padding=6)
        self.help_frame.pack(fill=tk.X, padx=5, pady=(6,2))

        self.help_rows = []
        help_lines = [
            ("🖱 左键拖拽峰圆点",  "水平移动峰位"),
            ("🖱 Shift+拖拽",      "上下拖动调整强度"),
            ("🖱 右键空白处",      "添加新峰"),
            ("🖱 右键峰上",        "删除该峰"),
            ("⌨ ← → 方向键",     "微调位置 (±0.01°)"),
            ("⌨ Shift+←→",       "微调强度 (±50)"),
            ("⌨ Ctrl+Z / Ctrl+Y","撤销 / 重做"),
            ("📝 双击表格",        "直接编辑数值"),
            ("📝 双击图谱空白",    "自动寻峰"),
            ("📂 拖放TXT到窗口",   "直接导入数据"),
        ]
        for icon_text, desc in help_lines:
            row = ttk.Frame(self.help_frame)
            row.pack(fill=tk.X, pady=1)
            ttk.Label(row, text=icon_text, font=('Microsoft YaHei', 8), foreground='#555',
                      width=18, anchor='w').pack(side=tk.LEFT)
            ttk.Label(row, text=desc, font=('Microsoft YaHei', 8, 'bold'),
                      foreground='#0066CC').pack(side=tk.LEFT)
            self.help_rows.append(row)

        # 折叠控制
        toggle_row = ttk.Frame(self.help_frame)
        toggle_row.pack(fill=tk.X, pady=(2,0))
        self.help_visible = tk.BooleanVar(value=True)
        ttk.Checkbutton(toggle_row, text="收起指南", variable=self.help_visible,
                         command=self._toggle_help).pack(side=tk.RIGHT)

        # 峰列表
        ttk.Label(right_frame, text="峰列表  (双击单元格编辑  |  单击选中峰)",
                  font=('Microsoft YaHei', 10, 'bold')).pack(pady=(4,2))

        cols = ('#', 'pos', 'h', 'fwhm', 'area')
        self.peak_tree = ttk.Treeview(right_frame, columns=cols, show='headings', height=10, selectmode='browse')
        self.peak_tree.heading('#', text='#');       self.peak_tree.column('#', width=25, anchor='center')
        self.peak_tree.heading('pos', text='位置(°)'); self.peak_tree.column('pos', width=72)
        self.peak_tree.heading('h', text='强度');      self.peak_tree.column('h', width=60)
        self.peak_tree.heading('fwhm', text='FWHM');  self.peak_tree.column('fwhm', width=65)
        self.peak_tree.heading('area', text='面积');    self.peak_tree.column('area', width=60)
        self.peak_tree.pack(fill=tk.BOTH, expand=True, padx=5)

        self.peak_tree.bind('<<TreeviewSelect>>', self._on_tree_select)
        self.peak_tree.bind('<Double-1>', self._on_tree_double)
        self.peak_tree.bind('<Delete>', lambda e: self._delete_selected())

        # 参数
        pf = ttk.LabelFrame(right_frame, text="参数", padding=5); pf.pack(fill=tk.X, padx=5, pady=4)

        ttk.Label(pf, text="峰形:").grid(row=0, column=0, sticky='w')
        ttk.Combobox(pf, textvariable=self.fit_method, width=14,
                     values=['Gaussian','Lorentzian','PseudoVoigt'], state='readonly').grid(row=0, column=1, sticky='w', padx=3)

        ttk.Label(pf, text="灵敏度:").grid(row=1, column=0, sticky='w', pady=(4,0))
        self.prom_var = tk.DoubleVar(value=3.0)
        ttk.Scale(pf, from_=0.5, to=20, variable=self.prom_var, orient=tk.HORIZONTAL,
                  command=lambda v: self.prom_lbl.configure(text=f"{float(v):.1f}%")).grid(row=1, column=1, sticky='ew', pady=(4,0))
        self.prom_lbl = ttk.Label(pf, text="3.0%"); self.prom_lbl.grid(row=1, column=2, padx=3, pady=(4,0))

        self.bg_var = tk.BooleanVar(); self.bg_var.trace_add('write', lambda *a: self._toggle_bg())
        ttk.Checkbutton(pf, text="背景扣除", variable=self.bg_var).grid(row=2, column=0, columnspan=2, sticky='w', pady=(4,0))
        self.residual_var = tk.BooleanVar()
        ttk.Checkbutton(pf, text="显示残差", variable=self.residual_var, command=self._toggle_residual).grid(row=3, column=0, columnspan=2, sticky='w')

        # 手动编辑区
        ef = ttk.LabelFrame(right_frame, text="选中峰参数  (可直接修改)", padding=5); ef.pack(fill=tk.X, padx=5, pady=4)

        self.edit_x = tk.StringVar(); self.edit_h = tk.StringVar(); self.edit_w = tk.StringVar()
        for j, (label, var) in enumerate([("位置:", self.edit_x), ("强度:", self.edit_h), ("FWHM:", self.edit_w)]):
            ttk.Label(ef, text=label).grid(row=j, column=0, sticky='w', pady=1)
            ttk.Entry(ef, textvariable=var, width=10).grid(row=j, column=1, sticky='w', padx=3, pady=1)
            var.trace_add('write', lambda *a, v=var, k=label: self._on_edit_field(v, k))
        ttk.Button(ef, text="应用修改", command=self._apply_edit).grid(row=3, column=0, columnspan=2, pady=(4,0))

        # 按钮
        bf = ttk.Frame(right_frame); bf.pack(fill=tk.X, padx=5, pady=6)
        for row_data in [
            [("📂 导入TXT", self.load_data), ("🔍 自动寻峰", self.auto_find_peaks)],
            [("📈 全谱拟合", self.fit_all), ("🗑 删除选中峰", self._delete_selected)],
            [("↩ 撤销", self.undo), ("↪ 重做", self.redo)],
            [("📋 导出峰参数(CSV)", self.export_csv), ("📄 导出精修数据(TXT)", self.export_txt)],
            [("🖼 导出图谱(PNG)", self.export_png), (None, None)],
        ]:
            r = ttk.Frame(bf); r.pack(fill=tk.X, pady=2)
            for label, cmd in row_data:
                if label is None: continue
                ttk.Button(r, text=label, command=cmd).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=1)

        # 状态栏
        self.status_var = tk.StringVar(value="就绪 — 请导入 XRD 数据文件 (.txt)")
        ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor='w').pack(fill=tk.X, side=tk.BOTTOM)

    # ── 拖放导入 ───────────────────────────
    def _setup_drag_drop(self):
        try:
            from tkinterdnd2 import DND_FILES  # 可选依赖
            self.root.drop_target_register(DND_FILES)
            self.root.dnd_bind('<<Drop>>', self._on_drop)
        except ImportError:
            pass  # tkinterdnd2 未安装，使用文件对话框即可

        # 使用 Windows shell 拖放 (仅 Windows)
        try:
            import ctypes
            self.root.drop_target_register('*')
        except:
            pass

    def _on_drop(self, event):
        path = event.data.strip('{}').strip()
        if os.path.isfile(path):
            self._do_load(path)

    # ── 加载数据 ───────────────────────────
    def load_data(self, filepath=None):
        if not filepath:
            filepath = filedialog.askopenfilename(filetypes=[("文本文件","*.txt *.dat *.csv"), ("所有","*.*")])
        if filepath: self._do_load(filepath)

    def _do_load(self, path):
        try:
            n = self.data.load(path)
            self.y_fit = None; self.r_factor = 0
            self.selected_peak = -1; self.detector.peaks.clear()
            self.undo_stack.clear(); self.redo_stack.clear(); self._push_undo()
            self._redraw(); self._update_peak_tree()
            self.file_label.configure(text=f"📄 {os.path.basename(path)}  ({n} 点)", foreground='#333')
            self.status_var.set(f"已加载: {os.path.basename(path)}  |  {n} 个数据点  |  点击「自动寻峰」开始")
        except Exception as e:
            messagebox.showerror("加载失败", str(e))

    # ── 事件 ────────────────────────────────
    def _connect_events(self):
        self.fig.canvas.mpl_connect('button_press_event', self._on_press)
        self.fig.canvas.mpl_connect('button_release_event', self._on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_motion)
        self.fig.canvas.mpl_connect('scroll_event', self._on_scroll)
        self.fig.canvas.mpl_connect('key_press_event', self._on_key)
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_hover)

    def _on_press(self, event):
        if event.inaxes != self.ax_main or self.data.x is None: return
        if event.dblclick: self._on_dbl(event); return
        if event.button == 3: self._on_right(event); return
        if event.button == 1:
            idx, _ = self.detector.find_nearest(event.xdata, x_tol=0.5)
            if idx >= 0:
                self.selected_peak = idx
                self.dragging_peak = idx
                self.drag_start = (event.xdata, event.ydata)
                # Shift = 调强度, 否则 = 调位置
                self.drag_mode = 'height' if (event.guiEvent and event.guiEvent.state & 0x1) else 'x'
                self._push_undo()
            else:
                self.selected_peak = -1; self.dragging_peak = -1
            self._update_peak_markers(); self._update_peak_tree(); self._update_edit_fields()

    def _on_release(self, event):
        if self.dragging_peak >= 0:
            p = self.detector.peaks[self.dragging_peak]
            self.status_var.set(f"峰 #{self.dragging_peak+1}: 位置={p['x']:.3f}°  强度={p['height']:.0f}")
        self.dragging_peak = -1; self.drag_mode = None; self.drag_start = None

    def _on_motion(self, event):
        if self.dragging_peak < 0 or event.inaxes != self.ax_main: return
        if event.xdata is None or event.ydata is None: return
        p = self.detector.peaks[self.dragging_peak]

        if self.drag_mode == 'x':
            dx = event.xdata - self.drag_start[0]
            p['x'] = np.clip(p['x'] + dx, self.data.x[0], self.data.x[-1])
        elif self.drag_mode == 'height':
            dy = event.ydata - self.drag_start[1]
            p['height'] = max(1, p['height'] + dy)
            p['area'] = p['height'] * p['fwhm'] * 1.064

        self.drag_start = (event.xdata, event.ydata)
        self._update_fit_preview(); self._update_fit_curve()
        self._update_peak_markers(); self._update_peak_tree(); self._update_edit_fields()

    def _on_right(self, event):
        if event.xdata is None: return
        self._push_undo()
        idx, _ = self.detector.find_nearest(event.xdata, x_tol=0.4)
        y_disp = self.data.y_corrected if self.bg_enabled else self.data.y_raw
        if idx >= 0:
            self.detector.remove_peak(idx)
            self.selected_peak = -1; self._update_fit_preview()
            self.status_var.set(f"已删除峰 #{idx+1}")
        else:
            rng = 0.3
            mask = (self.data.x >= event.xdata-rng) & (self.data.x <= event.xdata+rng)
            if mask.sum() > 2:
                h = y_disp[mask].max()
                x0 = self.data.x[mask][np.argmax(y_disp[mask])]
            else:
                x0 = event.xdata; h = y_disp[np.argmin(np.abs(self.data.x-x0))]
            self.detector.add_peak(x0, max(h, 1))
            self.status_var.set(f"已添加峰 @ {x0:.3f}°")
        self._redraw(); self._update_peak_tree()

    def _on_dbl(self, event):
        if event.xdata is None: return
        idx, _ = self.detector.find_nearest(event.xdata, x_tol=0.5)
        if idx < 0: self.auto_find_peaks()
        else:
            self.selected_peak = idx
            self._update_peak_markers(); self._update_peak_tree()
            self.peak_tree.selection_set(str(idx)); self.peak_tree.see(str(idx))

    def _on_hover(self, event):
        if event.inaxes != self.ax_main or self.data.x is None: return
        idx, d = self.detector.find_nearest(event.xdata or 0, x_tol=0.3)
        # 鼠标样式提示
        if self.dragging_peak >= 0: return

    def _on_scroll(self, event):
        if event.inaxes != self.ax_main or self.data.x is None: return
        s = 1.15 if event.button == 'up' else 0.85
        lo, hi = self.ax_main.get_ylim()
        mid = (lo+hi)/2; rng = (hi-lo)*s
        self.ax_main.set_ylim(mid-rng/2, mid+rng/2); self.canvas.draw_idle()

    def _on_key(self, event):
        if event.key == 'delete': self._delete_selected()
        elif event.key == 'ctrl+z': self.undo()
        elif event.key == 'ctrl+y': self.redo()
        elif event.key in ('left','right') and self.selected_peak >= 0:
            self._push_undo()
            step = 0.01 if event.key == 'right' else -0.01
            if event.guiEvent and event.guiEvent.state & 0x1:  # Shift
                self.detector.update_peak(self.selected_peak, height=max(1,
                    self.detector.peaks[self.selected_peak]['height'] + (50 if step>0 else -50)))
            else:
                self.detector.update_peak(self.selected_peak, x=np.clip(
                    self.detector.peaks[self.selected_peak]['x'] + step, self.data.x[0], self.data.x[-1]))
            self._update_peak_markers(); self._update_peak_tree(); self._update_edit_fields()
            self.canvas.draw()  # 立即刷新

    # ── 表格交互 ────────────────────────────
    def _on_tree_select(self, event):
        sel = self.peak_tree.selection()
        if sel:
            self.selected_peak = int(sel[0])
            self._update_peak_markers(); self._update_edit_fields()

    def _on_tree_double(self, event):
        """双击表格单元格 → 内联编辑"""
        col = self.peak_tree.identify_column(event.x)
        item = self.peak_tree.identify_row(event.y)
        if not col or not item: return
        col_idx = int(col.replace('#','')) - 1
        if col_idx < 0: return
        idx = int(item)
        p = self.detector.peaks[idx]
        keys = ['pos', 'h', 'fwhm']
        if col_idx-1 >= len(keys): return
        key = keys[col_idx-1]
        old_val = f"{p['x']:.3f}" if key=='pos' else f"{p['height']:.0f}" if key=='h' else f"{p['fwhm']:.3f}"

        bbox = self.peak_tree.bbox(item, col)
        self.edit_entry = ttk.Entry(self.peak_tree, width=10)
        self.edit_entry.place(x=bbox[0], y=bbox[1], width=bbox[2], height=bbox[3])
        self.edit_entry.insert(0, old_val); self.edit_entry.select_range(0, 'end'); self.edit_entry.focus()
        self.edit_entry.bind('<Return>', lambda e, i=idx, k=key: self._commit_edit(i, k))
        self.edit_entry.bind('<FocusOut>', lambda e: self._cancel_edit())
        self.edit_entry.bind('<Escape>', lambda e: self._cancel_edit())

    def _commit_edit(self, idx, key):
        if self.edit_entry is None: return
        try:
            val = float(self.edit_entry.get())
            self._push_undo()
            if key == 'pos':
                self.detector.update_peak(idx, x=np.clip(val, self.data.x[0], self.data.x[-1]))
            elif key == 'h':
                self.detector.update_peak(idx, height=max(1, val))
            elif key == 'fwhm':
                self.detector.update_peak(idx, fwhm=max(0.01, val))
            self._update_fit_preview()
        except ValueError:
            pass
        self._cancel_edit()
        self._redraw(); self._update_peak_tree()

    def _cancel_edit(self):
        if self.edit_entry:
            self.edit_entry.destroy(); self.edit_entry = None

    # ── 编辑字段 ────────────────────────────
    def _update_edit_fields(self):
        if self.selected_peak >= 0 and self.selected_peak < len(self.detector.peaks):
            p = self.detector.peaks[self.selected_peak]
            self.edit_x.set(f"{p['x']:.4f}")
            self.edit_h.set(f"{p['height']:.1f}")
            self.edit_w.set(f"{p['fwhm']:.4f}")
        else:
            self.edit_x.set(""); self.edit_h.set(""); self.edit_w.set("")

    def _on_edit_field(self, var, label):
        """编辑框内容变化时实时更新图谱"""
        if self.selected_peak < 0 or self.data.x is None: return
        try:
            val = float(var.get())
        except ValueError:
            return  # 还没输完数字, 跳过
        if label == "位置:":
            self.detector.update_peak(self.selected_peak,
                x=np.clip(val, self.data.x[0], self.data.x[-1]))
        elif label == "强度:":
            self.detector.update_peak(self.selected_peak, height=max(1, val))
        elif label == "FWHM:":
            self.detector.update_peak(self.selected_peak, fwhm=max(0.01, val))
        self._update_fit_preview()
        self._update_peak_markers(); self._update_peak_tree()
        self._update_fit_curve()
        self.canvas.draw()  # 立即刷新

    def _update_fit_preview(self):
        """用当前峰参数快速生成预览曲线（不重新拟合）"""
        if self.data.x is None or not self.detector.peaks:
            self.y_fit = None; self.r_factor = 0; return
        x = self.data.x
        y_preview = np.zeros_like(x)
        for p in self.detector.peaks:
            sigma = p['fwhm'] / 2.355
            y_preview += p['height'] * np.exp(-((x - p['x'])**2) / (2*sigma**2))
        bg = self.data.y_bg if self.bg_enabled else np.zeros_like(x)
        self.y_fit = y_preview + bg
        self.r_factor = 0  # 预览不是真实拟合, 清除R因子

    def _update_fit_curve(self):
        """轻量更新拟合曲线（不重绘全部）"""
        if self.y_fit is None:
            if self.fit_line:
                try: self.fit_line.remove()
                except: pass
                self.fit_line = None
            return
        if self.fit_line is None:
            self.fit_line, = self.ax_main.plot(self.data.x, self.y_fit,
                color=COLORS['fit'], lw=1.5, label='预览曲线')
        else:
            self.fit_line.set_data(self.data.x, self.y_fit)

    def _apply_edit(self):
        if self.selected_peak < 0: return
        self._push_undo()
        try:
            x = float(self.edit_x.get())
            h = float(self.edit_h.get())
            w = float(self.edit_w.get())
            self.detector.update_peak(self.selected_peak,
                x=np.clip(x, self.data.x[0], self.data.x[-1]),
                height=max(1, h), fwhm=max(0.01, w))
            self._update_fit_preview()
            self._redraw(); self._update_peak_tree()
            self.status_var.set(f"峰 #{self.selected_peak+1} 已更新")
        except ValueError:
            self.status_var.set("请输入有效数字")

    # ── 绘制 ────────────────────────────────
    def _redraw(self):
        ax = self.ax_main; ax.clear()
        if self.data.x is None:
            ax.set_title("请导入 XRD 数据文件 (.txt)"); self.canvas.draw(); return

        y_disp = self.data.y_corrected if self.bg_enabled else self.data.y_raw
        self.raw_line = ax.plot(self.data.x, y_disp, color=COLORS['raw'], lw=0.8, label='XRD 数据')[0]
        if self.bg_enabled:
            self.bg_line = ax.plot(self.data.x, self.data.y_bg, color=COLORS['bg'], lw=0.8, ls='--', label='背景基线')[0]
        if self.y_fit is not None:
            self.fit_line = ax.plot(self.data.x, self.y_fit, color=COLORS['fit'], lw=1.5, label='拟合曲线')[0]
        if self.show_residual and self.y_fit is not None:
            res = y_disp - self.y_fit
            y0 = ax.get_ylim()[0] if ax.get_ylim()[0] > -1e9 else -50
            ax.plot(self.data.x, res + y0, color=COLORS['residual'], lw=0.5, alpha=0.5)

        self.peak_lines.clear(); self.peak_dots.clear(); self.peak_labels.clear()
        for i, p in enumerate(self.detector.peaks):
            sel = (i == self.selected_peak)
            c = COLORS['selected'] if sel else COLORS['peak_line']
            lw = 2.0 if sel else 0.8
            y_at = y_disp[np.argmin(np.abs(self.data.x - p['x']))]
            vline = ax.axvline(p['x'], color=c, lw=lw, alpha=0.85)
            self.peak_lines.append(vline)
            dot = ax.plot(p['x'], p['height'], 'o', color=COLORS['peak_dot'],
                          markersize=8 if sel else 5, picker=True, pickradius=10, zorder=5)[0]
            self.peak_dots.append(dot)
            lbl = ax.annotate(f"{p['x']:.2f}\n{p['height']:.0f}", (p['x'], p['height']),
                               textcoords="offset points", xytext=(0, 14), fontsize=7,
                               color=c, ha='center', alpha=0.9)
            self.peak_labels.append(lbl)

        self.height_line = None
        if self.selected_peak >= 0 and self.selected_peak < len(self.detector.peaks):
            h = self.detector.peaks[self.selected_peak]['height']
            self.height_line = ax.axhline(h, color=COLORS['height_line'], lw=1, ls=':', alpha=0.7)

        ax.set_xlabel("2θ (度)", fontsize=11)
        ax.set_ylabel("强度 (counts)", fontsize=11)
        ax.set_title(f"XRD — {os.path.basename(self.data.filepath) if self.data.filepath else 'Data'}", fontsize=12)
        ax.legend(loc='upper right', fontsize=8, framealpha=0.6)
        ax.set_xlim(self.data.x[0], self.data.x[-1])
        ax.set_ylim(-y_disp.max()*0.03, y_disp.max()*1.18)
        self.canvas.draw()

    def _update_peak_markers(self):
        if self.data.x is None: return
        y_disp = self.data.y_corrected if self.bg_enabled else self.data.y_raw
        for i, p in enumerate(self.detector.peaks):
            if i >= len(self.peak_lines): break
            sel = (i == self.selected_peak)
            c = COLORS['selected'] if sel else COLORS['peak_line']
            self.peak_lines[i].set_xdata([p['x']]*2); self.peak_lines[i].set_color(c)
            self.peak_lines[i].set_linewidth(2.0 if sel else 0.8)
            self.peak_dots[i].set_data([p['x']], [p['height']])
            self.peak_dots[i].set_color(c); self.peak_dots[i].set_markersize(8 if sel else 5)
            self.peak_labels[i].set_position((p['x'], p['height']))
            self.peak_labels[i].set_text(f"{p['x']:.2f}\n{p['height']:.0f}")
            self.peak_labels[i].set_color(c)
        # 高度参考线
        if self.height_line:
            if self.selected_peak >= 0 and self.selected_peak < len(self.detector.peaks):
                self.height_line.set_ydata([self.detector.peaks[self.selected_peak]['height']]*2)
        elif self.selected_peak >= 0:
            self._redraw()  # 需要创建高度线
            return
        self.canvas.draw_idle()

    def _update_peak_tree(self):
        for item in self.peak_tree.get_children(): self.peak_tree.delete(item)
        for i, p in enumerate(self.detector.peaks):
            tag = 'sel' if i == self.selected_peak else ''
            self.peak_tree.insert('', 'end', iid=str(i),
                values=(i+1, f"{p['x']:.3f}", f"{p['height']:.0f}", f"{p['fwhm']:.3f}", f"{p['area']:.0f}"), tags=(tag,))
        self.peak_tree.tag_configure('sel', background='#FFDDDD')
        self._update_status()

    def _update_status(self):
        r_s = f" | R²={self.r_factor:.4f}" if self.r_factor > 0 else ""
        bg_s = " | 背景ON" if self.bg_enabled else ""
        self.status_var.set(f"峰数:{len(self.detector.peaks)} | 数据:{len(self.data.x) if self.data.x is not None else 0}点{r_s}{bg_s}")

    # ── 撤销/重做 ───────────────────────────
    def _push_undo(self):
        self.undo_stack.append({
            'peaks': copy.deepcopy(self.detector.peaks),
            'sel': self.selected_peak, 'bg': self.bg_enabled,
            'y_fit': self.y_fit.copy() if self.y_fit is not None else None,
            'r': self.r_factor
        })
        if len(self.undo_stack) > 30: self.undo_stack.pop(0)
        self.redo_stack.clear()

    def undo(self):
        if len(self.undo_stack) <= 1: self.status_var.set("无可撤销操作"); return
        self.redo_stack.append(self.undo_stack.pop())
        s = self.undo_stack[-1]; self._restore(s)
        self.status_var.set(f"已撤销 ({len(self.undo_stack)-1} 步)")

    def redo(self):
        if not self.redo_stack: self.status_var.set("无可重做操作"); return
        s = self.redo_stack.pop(); self.undo_stack.append(s)
        self._restore(s); self.status_var.set("已重做")

    def _restore(self, s):
        self.detector.peaks = copy.deepcopy(s['peaks'])
        self.selected_peak = s['sel']; self.bg_enabled = s['bg']
        self.bg_var.set(self.bg_enabled)
        self.y_fit = s['y_fit'].copy() if s['y_fit'] else None
        self.r_factor = s['r']
        self._redraw(); self._update_peak_tree(); self._update_edit_fields()

    # ── 功能 ────────────────────────────────
    def auto_find_peaks(self):
        if self.data.x is None: messagebox.showinfo("提示","请先导入数据"); return
        self._push_undo()
        y_disp = self.data.y_corrected if self.bg_enabled else self.data.y_raw
        prom = y_disp.max() * self.prom_var.get() / 100
        dist = max(3, len(self.data.x)//300)
        self.detector.auto_detect(self.data.x, y_disp, prominence=prom, distance=dist)
        self.selected_peak = -1; self.y_fit = None; self.r_factor = 0
        self._redraw(); self._update_peak_tree()
        self.status_var.set(f"自动寻峰完成 → {len(self.detector.peaks)} 个峰 | 灵敏度={self.prom_var.get():.1f}%")

    def fit_all(self):
        if self.data.x is None or not self.detector.peaks:
            messagebox.showinfo("提示","请先导入数据并寻峰"); return
        self._push_undo()
        y = self.data.y_raw
        bg = self.data.y_bg if self.bg_enabled else None
        yt, yc, _, r = self.fitter.fit(self.data.x, y, self.detector.peaks, bg=bg, method=self.fit_method.get())
        if yt is not None:
            self.y_fit = yt; self.r_factor = r
            self._redraw(); self._update_peak_tree(); self._update_edit_fields()
            self.status_var.set(f"拟合完成 | R²={r:.4f} | {len(self.detector.peaks)} 峰 | {self.fit_method.get()}")
        else:
            self.status_var.set("拟合失败 | 尝试减少峰数或更换峰形")

    def _toggle_help(self):
        show = self.help_visible.get()
        for row in self.help_rows:
            if show:
                row.pack(fill=tk.X, pady=1)
            else:
                row.pack_forget()

    def _toggle_bg(self):
        if self.bg_var.get():
            self.data.subtract_background(); self.bg_enabled = True
        else:
            self.bg_enabled = False
        self._update_fit_preview()
        self._redraw(); self._update_peak_tree()

    def _toggle_residual(self):
        self.show_residual = self.residual_var.get(); self._redraw()

    def _delete_selected(self):
        if self.selected_peak < 0: return
        self._push_undo()
        self.detector.remove_peak(self.selected_peak)
        self.selected_peak = -1; self.y_fit = None
        self._redraw(); self._update_peak_tree(); self._update_edit_fields()
        self.status_var.set("峰已删除")

    def export_csv(self):
        if not self.detector.peaks: messagebox.showinfo("提示","无峰数据"); return
        p = filedialog.asksaveasfilename(defaultextension='.csv', filetypes=[("CSV","*.csv")], initialfile='xrd_peaks.csv')
        if not p: return
        with open(p, 'w', newline='', encoding='utf-8-sig') as f:
            w = csv.writer(f)
            w.writerow(['#','位置(°)','强度','FWHM','面积','峰形'])
            for i, pk in enumerate(self.detector.peaks):
                w.writerow([i+1, f"{pk['x']:.4f}", f"{pk['height']:.1f}", f"{pk['fwhm']:.4f}", f"{pk['area']:.1f}", pk['type']])
        self.status_var.set(f"峰参数已导出: {os.path.basename(p)}")

    def export_txt(self):
        """导出完整精修数据 TXT"""
        if self.data.x is None: messagebox.showinfo("提示","无数据"); return
        p = filedialog.asksaveasfilename(defaultextension='.txt', filetypes=[("TXT","*.txt")], initialfile='xrd_refined_data.txt')
        if not p: return

        y_raw = self.data.y_raw
        y_bg = self.data.y_bg if self.bg_enabled else np.zeros_like(self.data.x)
        y_corr = self.data.y_corrected if self.bg_enabled else y_raw
        y_fit = self.y_fit if self.y_fit is not None else np.full_like(self.data.x, np.nan)
        residual = y_corr - y_fit

        header = "2Theta(deg)\tRawIntensity\tBackground\tCorrected\tFitted\tResidual"
        data = np.column_stack([self.data.x, y_raw, y_bg, y_corr, y_fit, residual])
        np.savetxt(p, data, fmt='%.4f', delimiter='\t', header=header)
        self.status_var.set(f"精修数据已导出: {os.path.basename(p)}")

    def export_png(self):
        p = filedialog.asksaveasfilename(defaultextension='.png', filetypes=[("PNG","*.png"),("PDF","*.pdf")], initialfile='xrd_refined.png')
        if not p: return
        self.fig.savefig(p, dpi=200, bbox_inches='tight'); self.status_var.set(f"已导出: {os.path.basename(p)}")

    def _on_close(self):
        self.root.quit(); self.root.destroy()

# ============================================================
def main():
    root = tk.Tk()
    app = XrdRefinerApp(root)
    # 自动加载或命令行参数
    data_path = None
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        data_path = sys.argv[1]
    else:
        p = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'TP.txt')
        if os.path.exists(p): data_path = p
    if data_path:
        try:
            app._do_load(data_path)
        except Exception:
            pass
    root.mainloop()

if __name__ == '__main__':
    main()
