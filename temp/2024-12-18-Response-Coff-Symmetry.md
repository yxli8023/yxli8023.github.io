---
title: 对称性约束响应系数
tags:  transport Group-Theory
layout: article
license: true
toc: true
key: a20241218
pageview: true
cover: /assets/images/Mma/fft-2.png
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

整理一下学习对称性约束响应系数的笔记
{:.info}
<!--more-->


最近在学习新东西以及看文章的时候涉及到了对称性约束响应系数的分析，虽然之前也简单的了解过，但始终没有细致的去推导过对称操作到底是怎么作用在系统上并约束响应系数的，这里借着一篇文献的中的结果，学习并整理一下这方面的笔记。这里主要是参考这里[Nonlinear Superconducting Magnetoelectric Effect](https://arxiv.org/abs/2404.18616)这篇文文章的分析结果，具体的一些物理细节以及响应函数可以参考这篇文章，在后面就直接从响应的物理量和外场进行分析了。


# 一阶响应

首先是体系磁化与Cooper对动量之间，在线性阶满足

$$
\begin{equation}
\delta M_a^{(1)}=\alpha_{ab}q_b
\end{equation}
$$

其中的$\delta M_a^{(1)}$是体系的磁化，$q_b$是Cooper的动量(其实这里可以不用局限于这里是Cooper对，只要认为这是一个物理量，要分析它在对称性下面的变换关系即可)。假如体系具有时间反演对称性$\mathcal{T}$，那么在$\mathcal{T}$的操作下$\delta M_a^{(1)}$是要发生反号的

$$
\begin{equation}
\begin{aligned}
{\color{blue}\mathcal{T}\delta M_a^{(1)}\mathcal{T}^{-1}}&=-\delta M_a^{(1)}\\
{\color{blue}\color{blue}\mathcal{T}(\alpha_{ab}q_b)\mathcal{T}^{-1}}&=-\delta M_a^{(1)}\qquad \mathcal{T}q_b\mathcal{T}^{-1}=-q_b\\
\alpha_{ab}{\color{red}(-q)}&=-\delta M_a^{(1)}=-\alpha_{ab}q_b\Longrightarrow \alpha_{ab}\neq 0 
\end{aligned}
\end{equation}
$$

因此可以发现当系统具有时间反演对称性的时候，一阶响应系数非零，即不会被时间反演对称性禁止。

接下来再看空间反演对称性$\mathcal{P}$，首先要知道的是在空间反演操作下动量$q_b$会反号，而磁化是不会反号的，从而有

$$
\begin{equation}
\begin{aligned}
{\color{blue}\mathcal{P}\delta M_a^{(1)}\mathcal{P}^{-1}}&=\delta M_a^{(1)}\\
{\color{blue}\color{blue}\mathcal{P}(\alpha_{ab}q_b)\mathcal{P}^{-1}}&=\delta M_a^{(1)}\qquad \mathcal{P}q_b\mathcal{P}^{-1}=-q_b\\
\alpha_{ab}{\color{red}(-q)}&=\delta M_a^{(1)}=\alpha_{ab}q_b\Longrightarrow {\color{teal}\alpha_{ab}\equiv 0 }
\end{aligned}
\end{equation}
$$

因此如果体系除了时间反演对称性，还存在空间反演对称性$\mathcal{P}$，此时是不会存在一阶响应的。

# 二阶响应

接下来分析二阶响应

$$
\begin{equation}
\delta M_c^{(2)}=\chi^{c}_{ab}q_aq_b
\end{equation}
$$

首先还是分析时间反演对称性$\mathcal{T}$

$$
\begin{equation}
\begin{aligned}
{\color{blue}\mathcal{T}\delta M_c^{(2)}\mathcal{T}^{-1}}&=-\delta M_a^{(2)}\\
{\color{blue}\color{blue}\mathcal{T}(\chi^c_{ab}q_aq_b)\mathcal{T}^{-1}}&=-\delta M_a^{(2)}\qquad \mathcal{T}q_i\mathcal{T}^{-1}=-q_i\\
\chi^c_{ab}{\color{red}(-q_a)(-q_b)}&=-\delta M_a^{(2)}=-\chi^c_{ab}q_aq_b\Longrightarrow {\color{teal}\chi^c_{ab}\equiv 0} 
\end{aligned}
\end{equation}
$$

可以发现在具有时间反演对称性时二阶响应不存在$\chi^c_{ab}\equiv 0$。

而对于空间反演对称性$\mathcal{P}$则有

$$
\begin{equation}
\begin{aligned}
{\color{blue}\mathcal{P}\delta M_c^{(2)}\mathcal{P}^{-1}}&=\delta M_a^{(2)}\\
{\color{blue}\color{blue}\mathcal{P}(\chi^c_{ab}q_aq_b)\mathcal{P}^{-1}}&=\delta M_a^{(2)}\qquad \mathcal{P}q_i\mathcal{P}^{-1}=-q_i\\
\chi^c_{ab}{\color{red}(-q_a)(-q_b)}&=\delta M_a^{(2)}=\chi^c_{ab}q_aq_b\Longrightarrow {\color{teal}\chi^c_{ab}\neq 0} 
\end{aligned}
\end{equation}
$$

因此在空间反演对称性下，体系是具有二阶响应的，也就是如果该体系破坏了时间反演对称性$\mathcal{T}$，但具有空间反演对称性$\mathcal{P}$，那么体系则只会存在二阶响应。

综合上面的分析，可以发现如果体系同时具有$\mathcal{PT}$联合操作，那么体系的一阶响应$\alpha_{ab}\equiv0$以及二阶响应$\chi^c_{ab}\equiv$都是被对称性禁戒掉的。因此可以得到下面的一个表


<div align="center">

|       | $\mathcal{P}$   | $\mathcal{T}$   | $\mathcal{P} \mathcal{T}$ |
|-------|-----------------|-----------------|---------------------------|
| $M_x$ | ✔(2)            | ✔(1)            | ✗                         |
| $M_y$ | ✔(2)            | ✔(1)            | ✗                         |
| $M_z$ | ✔(2)            | ✔(1)            | ✗                         |

</div>


# 晶体对称性约束



# 参考文献

- [Nonlinear Superconducting Magnetoelectric Effect](https://arxiv.org/abs/2404.18616)




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

