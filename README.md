# XRD 数据精修工具

交互式 XRD 图谱寻峰、拟合、精修工具。**桌面版 EXE 免安装** + **Web 在线版**。

---

## 快速开始

### 🖥 桌面版（推荐）
> 下载即用，无需安装 Python

**[⬇ 下载 EXE](https://github.com/Chile-Dragon/xrd-refiner/releases)** → 双击 `XRD_Refiner.exe`

### 🌐 Web 在线版
> 浏览器直接使用

🔗 [Streamlit Cloud](https://xrd-refiner.streamlit.app/)（部署后生效）

本地启动：
```bash
pip install -r requirements.txt
streamlit run xrd_web.py
```

---

## 功能

| 功能 | 桌面版 | Web 版 |
|------|:---:|:---:|
| 上传 TXT/CSV/DAT 数据 | ✅ | ✅ |
| 自动寻峰（可调灵敏度） | ✅ | ✅ |
| 拖拽峰位 / 调整强度 | ✅ 鼠标拖拽 | ✅ 编辑面板 |
| 手动添加/删除峰 | ✅ | ✅ |
| Gaussian / Lorentzian / PseudoVoigt 拟合 | ✅ | ✅ |
| 背景扣除 / 残差显示 | ✅ | ✅ |
| 实时预览曲线 | ✅ | ✅ |
| 撤销 / 重做 | ✅ | ✅ |
| 导出峰参数 CSV | ✅ | ✅ |
| 导出精修数据 TXT | ✅ | ✅ |
| 导出图谱 PNG | ✅ | ✅ |
| 单峰组分展示 | — | ✅ |
| 快捷键微调 | ✅ 方向键 | ✅ 按钮 |

---

## 桌面版使用说明

### 操作指南

| 鼠标 | 操作 |
|------|------|
| 🖱 **左键拖拽峰圆点** | 水平移动峰位 |
| 🖱 **Shift + 拖拽** | 上下拖动调整强度 |
| 🖱 **右键空白处** | 添加新峰 |
| 🖱 **右键峰标记** | 删除该峰 |
| 📝 **双击图谱空白** | 自动寻峰 |
| 📝 **双击表格** | 直接编辑参数 |

| 键盘 | 操作 |
|------|------|
| ← → | 微调峰位 (±0.01°) |
| Shift + ← → | 微调强度 (±50) |
| Delete | 删除选中峰 |
| Ctrl+Z / Ctrl+Y | 撤销 / 重做 |

---

## 数据格式

两列数据（2θ 角度 + 强度），Tab/逗号/空格分隔均可，自动跳过标题行。

```
5.0000  1792.0000
5.0123  1746.0000
5.0245  1707.0000
...
```

---

## 环境要求

| 依赖 | 版本 |
|------|------|
| Python | 3.9+ |
| numpy | ≥1.24 |
| scipy | ≥1.10 |
| matplotlib | ≥3.5 |
| streamlit | ≥1.28 (仅 Web 版) |
| plotly | ≥5.14 (仅 Web 版) |

### 安装

```bash
# 一键安装
pip install -r requirements.txt

# 仅桌面版
pip install numpy scipy matplotlib
```

### 常见问题

| 问题 | 解决 |
|------|------|
| `ModuleNotFoundError` | `pip install -r requirements.txt` |
| Web 版无法导出 PNG | `pip install kaleido` |
| macOS tkinter 报错 | `brew install python-tk@3.11` |

---

## 打包为 EXE

```bash
pip install pyinstaller
pyinstaller xrd_refiner.spec --clean --noconfirm
```
