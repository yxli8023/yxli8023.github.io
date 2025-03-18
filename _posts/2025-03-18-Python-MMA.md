---
title: 利用Python实现Mathematica的配色方案
tags:  Code Python Mathematica
layout: article
license: true
toc: true
key: a20250318
pageview: true
cover: /assets/images/python/Chern-couple.png
header:
  theme: dark
  background: 'linear-gradient(to right, rgb(101, 78, 163), rgb(234, 175, 200))'
article_header:
  type: overlay
  theme: dark
  background_color: false
  background_image: 
    gradient: 'linear-gradient(to right, rgb(101, 78, 163), rgb(234, 175, 200))'
    image: false
aside:
    toc: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
该笔记主要整理了如何在Python中实现Mathematica的配色
{:.info}
<!--more-->
# 前言
不得不说Mathematica中的一些颜色还是很好看的，但是当数据文件较大的时候用Mathematica绘图就会有点卡顿，不如用Python来的方便一些，这里就整理一下如何将Mathematica中的配色数据导出来，然后再加载到Python中用来绘图

- 首先将你喜欢的颜色数据导出来，这里可以通过改变变量dk来控制颜色间隔的大小
```fortran
  dk = 1/255;
  colorFunction = ColorData["LightTemperatureMap", "ColorFunction"];
  colorValues = Table[colorFunction[x], {x, 0, 1, dk}];
  Export["LightTemperatureMap.dat", colorValues, "Table"];
```

- 下面就将这个颜色数据导入到Python中创建自己的colorbar
```python
with open("LightTemperatureMap.dat", "r") as file:
        lines = file.readlines()
    # 解析颜色数据
    colors = []
    for line in lines:
        # 去掉多余的字符（例如 "RGBColor[" 和 "]"）
        line = line.replace("RGBColor[", "").replace("]", "").strip()
        # 按逗号分割，提取 R, G, B 值
        r, g, b = map(float, line.split(","))
        colors.append((r, g, b))

    # 将颜色数据转换为 matplotlib 的 colormap 格式
    cmap = LinearSegmentedColormap.from_list("LightTemperatureMap", colors)
```

- 最后就在绘制密度图的时候直接使用这个颜色配置即可

# 完整代码
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
plt.rc('font', family='Times New Roman')
config = {
"font.size": 35,
"mathtext.fontset":'stix',
"font.serif": ['SimSun'],
}
rcParams.update(config) # Latex 字体设置
#-----------------------------------------------------------------------------------------------
def plot_Chern_Number_2D_couple():
    dataname1 = "Phase-V0-delta0-HP.dat"    # 该数据有三列,分别是[x0,y0,z0]
    dataname2 = "Phase-V0-delta0-HM.dat"
    picname = "Chern-couple.png"
    da = np.loadtxt(dataname1) 
    da2 = np.loadtxt(dataname2) 
    x0 = da[:, 0]  # delta0 数据
    y0 = da[:, 1]  # d0 数据
    z0 = np.array(da[:, 2])  # Chern 数数据
    z1 = np.array(da2[:, 2])  # Chern 数数据

    # 获取 delta0 和 d0 的范围
    delta0_min, delta0_max = np.min(x0), np.max(x0)
    d0_min, d0_max = np.min(y0), np.max(y0)

    # 将 z0 数据重塑为二维数组
    xn = int(np.sqrt(len(x0)))
    z0 = z0.reshape(xn, xn) + z1.reshape(xn, xn)

    # 创建绘图
    plt.figure(figsize=(10, 10))
    with open("LightTemperatureMap.dat", "r") as file:
        lines = file.readlines()
    # 解析颜色数据
    colors = []
    for line in lines:
        # 去掉多余的字符（例如 "RGBColor[" 和 "]"）
        line = line.replace("RGBColor[", "").replace("]", "").strip()
        # 按逗号分割，提取 R, G, B 值
        r, g, b = map(float, line.split(","))
        colors.append((r, g, b))

    # 将颜色数据转换为 matplotlib 的 colormap 格式
    cmap = LinearSegmentedColormap.from_list("LightTemperatureMap", colors)
    # sc = plt.imshow(z0, interpolation = 'bilinear', cmap="RdYlBu", origin='lower', extent=[d0_min, d0_max, delta0_min, delta0_max])
    # sc = plt.imshow(z0, cmap="RdYlBu", origin='lower', extent=[d0_min, d0_max, delta0_min, delta0_max])
    sc = plt.imshow(z0, interpolation = 'bilinear',cmap = cmap, origin='lower', extent=[d0_min, d0_max, delta0_min, delta0_max])

    # 添加 colorbar
    cb = plt.colorbar(sc, fraction = 0.04, ticks = [np.min(z0), np.max(z0)], extend='both')
    cb.ax.set_title(r"$C$", fontsize=30)
    cb.ax.tick_params(size=1)
    cb.ax.set_yticklabels([format(np.min(z0), ".1f"), format(np.max(z0), ".1f")])

    # 设置坐标轴标签
    font2 = {'family': 'Times New Roman', 'weight': 'normal', 'size': 40}
    plt.xlabel(r"$\Delta_0$", font2)
    plt.ylabel(r"$V_0$", font2)
    plt.title(r"$H^++H^-$")

    # 设置坐标轴刻度
    plt.yticks(np.linspace(delta0_min, delta0_max, 5), fontproperties = 'Times New Roman', size=30)  # 显示 5 个刻度
    plt.xticks(np.linspace(d0_min, d0_max, 5), fontproperties = 'Times New Roman', size=30)  # 显示 5 个刻度
    plt.tick_params(axis = 'x', width=1, length=8, direction = "in")
    plt.tick_params(axis = 'y', width=1, length=8, direction = "in")

    # 设置横纵方向比例一致
    ax = plt.gca()
    ax.set_aspect('auto')  # 自动调整比例
    # 或者手动设置比例
    # ax.set_aspect((delta0_max - delta0_min) / (d0_max - d0_min))  # 根据数据范围设置比例

    # 设置坐标轴样式
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["right"].set_linewidth(1.5)
    ax.spines["top"].set_linewidth(1.5)

    # 保存图像
    plt.savefig(picname, dpi=300, bbox_inches = 'tight')
    plt.close()
```

![png](/assets/images/python/Chern-couple.png)


# 文件下载
所有使用到的数据和绘图程序[点击这里下载](/assets/data/AA.zip)

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

<table>
  <tr>
    <!-- 图片单元格 -->
    <td style="width: 300px; height: 300px; text-align: center; vertical-align: middle; border: 1px solid #ccc; border-radius: 8px;">
      <img src="/assets/images/qrcode.jpg" alt="QR Code" width="300px" height="300px" style="border-radius: 8px;">
    </td>
    <!-- 文字单元格 -->
    <td style="width: 300px; height: 300px; text-align: center; vertical-align: middle; padding-left: 20px; border: 1px solid #ccc; border-radius: 8px;">
      <div>
        <h4 style="margin: 0;">Email</h4>
        <p style="margin: 5px 0;">yxli406@gmail.com</p>
      </div>
    </td>
  </tr>
</table>

