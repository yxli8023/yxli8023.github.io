---
title: 转角石墨烯摩尔纹
tags: Code Mathematica
layout: article
license: true
toc: true
key: a20210208a
pageview: true
cover: /assets/images/Mma/moire-1.png
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
aside:
    toc: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
从转角石墨烯发现超导以来,转角已经成为了一个非常火热的研究方向,这里就简单的利用Mathematica实现转角石墨烯中的摩尔纹.
<!--more-->
作为凝聚态物理的研究生,转角石墨烯的火热让人感受颇深,这里我就想利用程序来看看转角石墨烯中的摩尔纹是怎么随着转角而变化的.
# code
这里直接上代码和结果,这里第一种方案是利用点位置来作图

![png](/assets/images/Mma/moire-1.png)

![png](/assets/images/Mma/moire-1-1.png)

第二种方案是直接使用石墨烯的格点坐标来绘制

![png](/assets/images/Mma/moire-2.png)

![png](/assets/images/Mma/moire-2-1.png)

# 四方点阵转角
既然转角很火,那么四方点阵也自然可以研究转角,比如最近在转角双层铜基超导上发现拓扑的存在

![png](/assets/images/Mma/moire-3.png)

# 代码下载
上面的代码可以[点击这里下载](/assets/data/2021-02-08.nb)

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