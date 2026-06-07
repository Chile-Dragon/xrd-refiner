# XRD 数据精修工具

交互式 XRD 图谱寻峰、拟合、精修工具，支持桌面版和 Web 版。

---

## 环境要求

| 依赖 | 最低版本 | 说明 |
|------|---------|------|
| Python | 3.9+ | [python.org](https://www.python.org/downloads/) 下载 |
| pip | 21.0+ | 随 Python 自动安装 |

### 依赖库

| 库 | 用途 | 版本要求 |
|----|------|---------|
| numpy | 数值计算 | ≥1.24 |
| scipy | 寻峰/拟合算法 | ≥1.10 |
| matplotlib | 桌面版绘图 | ≥3.5 |
| streamlit | Web 框架 | ≥1.28 |
| plotly | Web 交互图表 | ≥5.14 |
| kaleido | PNG 导出（可选） | ≥0.2.1 |

### 一键安装

```bash
# Web 版（含桌面版依赖）
pip install -r requirements.txt

# 仅桌面版
pip install numpy scipy matplotlib

# PNG 导出支持（可选）
pip install kaleido
```

---

## 安装与运行

### Windows

```bash
# 1. 安装 Python 3.9+（勾选 "Add Python to PATH"）
# 2. 打开命令提示符 (cmd) 或 PowerShell
# 3. 克隆项目
git clone https://github.com/Chile-Dragon/xrd-refiner.git
cd xrd-refiner

# 4. 安装依赖
pip install -r requirements.txt

# 5. 启动
streamlit run xrd_web.py          # Web 版
python xrd_refiner.py             # 桌面版
```

### macOS

```bash
# 安装 Homebrew（如果没有）
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# 安装 Python
brew install python@3.11

# 克隆并安装
git clone https://github.com/Chile-Dragon/xrd-refiner.git
cd xrd-refiner
pip3 install -r requirements.txt

# 启动
streamlit run xrd_web.py
```

### Linux (Ubuntu/Debian)

```bash
# 安装 Python 和 pip
sudo apt update
sudo apt install python3 python3-pip python3-tk git -y

# 克隆项目
git clone https://github.com/Chile-Dragon/xrd-refiner.git
cd xrd-refiner

# 安装依赖
pip3 install -r requirements.txt

# 启动
streamlit run xrd_web.py
```

---

## Web 版（在线使用）

🔗 公网地址: [Streamlit Cloud](https://xrd-refiner.streamlit.app/)（部署后生效）

### 本地启动
```bash
streamlit run xrd_web.py
# 浏览器自动打开 http://localhost:8501
```

### 功能
- 📂 上传 TXT/CSV/DAT 格式 XRD 数据
- 🔍 自动寻峰（可调灵敏度）
- ✏️ 点击峰 → 编辑位置/强度/FWHM → 实时预览
- 📈 三种峰形拟合 (Gaussian / Lorentzian / PseudoVoigt)
- 📊 单峰组分展示 + 残差分析
- 📥 导出峰参数 CSV + 精修数据 TXT + 图谱 PNG
- ↩ 撤销/重做

---

## 桌面版

```bash
python xrd_refiner.py
```

支持鼠标拖拽峰位、Shift+拖拽调整强度、方向键微调、双击表格编辑。

---

## 数据格式

支持两列数据（2θ角度 + 强度），Tab/逗号/空格分隔均可，自动跳过标题行。

示例：
```
5.0000  1792.0000
5.0123  1746.0000
5.0245  1707.0000
...
```

---

## 常见问题

**Q: 启动报错 `ModuleNotFoundError`？**
```bash
pip install -r requirements.txt
```

**Q: Web 版无法导出 PNG？**
```bash
pip install kaleido
```

**Q: macOS 桌面版报错 tkinter？**
```bash
brew install python-tk@3.11
```

**Q: 打包为 exe？**
```bash
pip install pyinstaller
pyinstaller --onefile --windowed xrd_refiner.py
```
