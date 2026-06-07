# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller 精简打包配置 — xrd_refiner.py
策略: 排除未使用的 matplotlib后端 / scipy子模块 / tkinter测试
"""

block_cipher = None

# ---- 要排除的模块 ----
EXCLUDES = [
    # matplotlib 不需要的后端（已验证安全）
    'matplotlib.backends.backend_qt5agg',
    'matplotlib.backends.backend_qt5',
    'matplotlib.backends.backend_qtagg',
    'matplotlib.backends.backend_gtk3agg',
    'matplotlib.backends.backend_gtk3',
    'matplotlib.backends.backend_gtk4agg',
    'matplotlib.backends.backend_gtk4',
    'matplotlib.backends.backend_wx',
    'matplotlib.backends.backend_wxagg',
    'matplotlib.backends.backend_webagg',
    'matplotlib.backends.backend_webagg_core',
    'matplotlib.backends.backend_nbagg',
    'matplotlib.backends.backend_macosx',
    'matplotlib.backends.backend_cairo',
    'matplotlib.backends.backend_template',
    'matplotlib.backends.backend_pgf',
    'matplotlib.backends.backend_ps',
    'matplotlib.backends.backend_svg',
    'matplotlib.backends.backend_pdf',
]

a = Analysis(
    ['xrd_refiner.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=[
        'scipy.signal._peak_finding',
        'scipy.signal._peak_finding_utils',
        'scipy.optimize._minpack',
        'scipy.optimize._lbfgsb',
        'scipy.special',
        'scipy.special._ufuncs',
        'numpy.core._methods',
        'numpy.lib.format',
        'numpy.linalg',
        'matplotlib.backends.backend_tkagg',
        'tkinter.filedialog',
        'tkinter.ttk',
        'secrets',
        'tkinter.messagebox',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=EXCLUDES,
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='XRD_Refiner',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,       # Windows 无 strip
    upx=False,          # UPX 未安装
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,      # 正式版: 无控制台
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,
)
