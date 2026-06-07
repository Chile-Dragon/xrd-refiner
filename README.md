# XRD 数据精修工具

交互式 XRD 图谱寻峰、拟合、精修工具，支持桌面版和 Web 版。

## Web 版（在线使用）

🔗 **一键打开**: [Streamlit Cloud](https://xrd-refiner.streamlit.app/)（部署后生效）

### 功能
- 📂 上传 TXT/CSV/DAT 格式 XRD 数据
- 🔍 自动寻峰（可调灵敏度）
- ✏️ 点击峰 → 编辑位置/强度/FWHM → 实时预览
- 📈 全谱高斯多峰拟合，显示 R²
- 📥 导出峰参数 CSV + 精修数据 TXT
- ↩ 撤销/重做

### 本地运行
```bash
pip install -r requirements.txt
streamlit run xrd_web.py
```

## 桌面版

```bash
pip install numpy scipy matplotlib
python xrd_refiner.py
```

支持鼠标拖拽峰位、Shift+拖拽调整强度、方向键微调、双击表格编辑。

## 数据格式

支持两列数据（2θ角度 + 强度），Tab/逗号/空格分隔均可，自动跳过标题行。
