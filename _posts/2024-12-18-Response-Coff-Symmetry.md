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

综合上面的分析，可以发现如果体系同时具有$\mathcal{PT}$联合操作，那么体系的一阶响应$\alpha_{ab}\equiv0$以及二阶响应$\chi^c_{ab}\equiv0$都是被对称性禁戒掉的。因此可以得到下面的一个表


<div align="center">

|       | $\mathcal{P}$   | $\mathcal{T}$   | $\mathcal{P} \mathcal{T}$ |
|-------|-----------------|-----------------|---------------------------|
| $M_x$ | ✔(2)            | ✔(1)            | ✗                         |
| $M_y$ | ✔(2)            | ✔(1)            | ✗                         |
| $M_z$ | ✔(2)            | ✔(1)            | ✗                         |

</div>


# 晶体对称性约束

## 一阶响应

前面提到的都是比较简单的对称性，这里考虑晶体对称性对响应系数的约束。假设系统具有绕着$z$轴的$\mathcal{C}_4$对称操作，这是一个二维平面上逆时针旋转$90^\degree$的操作，对动量$\boldsymbol{q}=(q_x,q_y)$的操作为

$$
\begin{equation}
\mathcal{C}_4 : \mathbf{q} \to \mathbf{q}' = R \mathbf{q}, \quad R = 
\begin{pmatrix}
0 & -1 \\
1 & 0
\end{pmatrix}.
\end{equation}
$$

相对应的可以得到

$$
\begin{equation}
q_x' = -q_y, \quad q_y' = q_x.
 \end{equation}
$$

对于磁化 $\delta M_a$，它在$\mathcal{C}_4$操作下的变换规则与动量相同，具体为：

$$
\begin{equation}
\mathcal{C}_4 : \delta M_a \to \delta M_a' = R \delta M_a.
 \end{equation}
$$


响应关系由以下方程描述：

$$
\begin{equation}
\delta M_a^{(1)} = \alpha_{ab} q_b,
 \end{equation}
$$

其中 $\alpha_{ab}$ 是响应系数张量(因为这里先考虑一阶响应，所以这里就是一个矩阵)。

在$\mathcal{C}_4$旋转对称操作下：

- $\delta M_a$ 和 $q_b$ 分别按照旋转矩阵 $R$ 变换。
- 变换后的响应关系可以写为：

$$
\begin{equation}
\delta M_a'^{(1)} = \alpha_{ab}' q_b',
 \end{equation}
$$

将$\mathcal{C}_4$操作作用到原始响应关系上，磁化和动量的变换可以代入：

$$
\begin{equation}
R_{ac} {\color{blue}\delta M_c^{(1)}} = \alpha_{ab}' R_{bd} q_d.
 \end{equation}
$$

将原始的响应关系 ${\color{blue}\delta M_c^{(1)} = \alpha_{ce} q_e}$ 代入上式：

$$
\begin{equation}
R_{ac} {\color{blue}\alpha_{ce} q_e} = \alpha_{ab}' R_{bd} q_d.
 \end{equation}
$$

由于这个方程必须对任意动量 $q$ 成立，我们可以去掉动量 $q$，得到：

$$
\begin{equation}
R_{ac} \alpha_{ce} = \alpha_{ab}' R_{be}.
 \end{equation}
$$

这说明响应系数 $\alpha_{ab}$ 必须满足如下关系：

$$
\begin{equation}
{\color{blue}\alpha_{ab}' = R_{ac} \alpha_{cd} R_{bd}^{-1}.}
 \end{equation}
$$

对于$\mathcal{C}_4$操作，旋转矩阵 $R$ 的逆等于它的转置，即 ${\color{blue}R^{-1} = R^T}$。因此，上式可以写为：

$$
\begin{equation}
\alpha_{ab}' = R_{ac} \alpha_{cd} {\color{blue}R_{db}}.
 \end{equation}
$$
{:.success}

现在分析响应张量 $\alpha_{ab}$ 的对称性。具体来说：
- 旋转矩阵 $R$ 在$\mathcal{C}_4$操作下为：

$$
\begin{equation}
R = 
\begin{pmatrix}
0 & -1 \\
1 & 0
\end{pmatrix}.
 \end{equation}
$$

将 $R$ 代入上述方程：

$$
\begin{equation}
\alpha_{ab} = R_{ac} \alpha_{cd} R_{bd}.
 \end{equation}
$$

展开矩阵乘法，逐项计算，得到：

$$
\begin{equation}
\begin{pmatrix}
\alpha_{xx} & \alpha_{xy} \\
\alpha_{yx} & \alpha_{yy}
\end{pmatrix}
=
\begin{pmatrix}
0 & -1 \\
1 & 0
\end{pmatrix}
\begin{pmatrix}
\alpha_{xx} & \alpha_{xy} \\
\alpha_{yx} & \alpha_{yy}
\end{pmatrix}
\begin{pmatrix}
0 & 1 \\
-1 & 0
\end{pmatrix}.
 \end{equation}
$$

我们逐步计算右侧矩阵乘法：

1. **第一步：中间乘法**：

$$
\begin{equation}
\begin{pmatrix}
0 & -1 \\
1 & 0
\end{pmatrix}
\begin{pmatrix}
\alpha_{xx} & \alpha_{xy} \\
\alpha_{yx} & \alpha_{yy}
\end{pmatrix}
=
\begin{pmatrix}
-\alpha_{yx} & -\alpha_{yy} \\
\alpha_{xx} & \alpha_{xy}
\end{pmatrix}.
 \end{equation}
$$

2. **第二步：与 $R^T$ 相乘**：

$$
\begin{equation}
\begin{pmatrix}
-\alpha_{yx} & -\alpha_{yy} \\
\alpha_{xx} & \alpha_{xy}
\end{pmatrix}
\begin{pmatrix}
0 & 1 \\
-1 & 0
\end{pmatrix}
=
\begin{pmatrix}
\alpha_{yy} & -\alpha_{yx} \\
-\alpha_{xy} & \alpha_{xx}
\end{pmatrix}.
 \end{equation}
$$

因此，$\alpha_{ab}$ 必须满足：

$$
\begin{equation}
\alpha_{xx} = \alpha_{yy}, \quad \alpha_{xy} = -\alpha_{yx}.
\end{equation}
$$

---


因此在$\mathcal{C}_4$旋转对称性下，响应张量 $\alpha_{ab}$ 的形式被约束为：

$$
\begin{equation}
\alpha_{ab} =
\begin{pmatrix}
\alpha & \beta \\
-\beta & \alpha
\end{pmatrix},
 \end{equation}
$$
其中：
- $\alpha$ 是对角元，表示各向同性的部分。
- $\beta$ 是反对称元，表示与方向相关的部分。

这种形式意味着：
- $\alpha_{xx} = \alpha_{yy} = \alpha$。
- $\alpha_{xy} = -\alpha_{yx} = \beta$。

这与$\mathcal{C}_4$对称性要求的旋转不变性相一致。


- **对角元 $\alpha$**：表示各向同性响应，即沿 $x$ 和 $y$ 方向的响应是相等的。
- **非对角元 $\beta$**：表示反对称响应，即系统具有某种与方向相关的奇性，比如霍尔效应类的现象。

同样可以拿代码展示一下这个约束过程

![png](/assets/images/GroupTheory/sym-1.png)


## 二阶响应

下面再来考虑晶体对称操作对二阶响应的约束，这里考虑的二阶响应关系为

$$
\begin{equation}
\delta M_c^{(2)} = \chi_{ab}^c q_a q_b,
 \end{equation}
$$

其中：
- $\chi_{ab}^c$ 是三阶张量，表示系统的二阶响应系数
- $q_a q_b$ 是动量的二阶乘积项


在$\mathcal{C}_4$变换下，二阶响应关系 $\delta M_c^{(2)} = \chi_{ab}^c q_a q_b$ 变为：

$$
\begin{equation}
\delta M_c'^{(2)} = \chi_{ab}'^c {\color{red}q_a'}{\color{blue} q_b'},
 \end{equation}
$$

将动量变换 ${\color{red}q_a' = R_{ae} q_e},{\color{green}q_b' = R_{bf} q_f}$ 和磁化变换 ${\color{blue}\delta M_c' = R_{cd} \delta M_d}$ 代入原始关系：

$$
\begin{equation}
R_{cd}  {\color{blue}\delta M_d^{(2)}} = \chi_{ab}'^c {\color{red}(R_{ae} q_e)}{\color{green}(R_{bf} q_f)}.
 \end{equation}
$$


由原始关系 ${\color{blue}\delta M_d^{(2)} = \chi_{ef}^d q_e q_f}$，代入上式可得：

$$
\begin{equation}
R_{cd} {\color{blue}\chi_{ef}^d {q_e q_f} }= \chi_{ab}'^c R_{ae} R_{bf} q_e q_f.
 \end{equation}
$$

因为这必须对任意 $q_e q_f$ 成立，我们可以去掉 $q$ 并比较系数

$$
\begin{equation}
R_{cd} {\color{black}\chi_{ef}^d \sout{q_e q_f} }= \chi_{ab}'^c R_{ae} R_{bf} \sout{q_e q_f}.
 \end{equation}
$$

得到：

$$
\begin{equation}
R_{cd}\chi_{ef}^d  = \chi_{ab}'^cR_{ae} R_{bf} .
 \end{equation}
$$

将等式右边乘以$R_{bf}^{-1}R_{ae}^{-1}$则得到

$$
\begin{equation}
R_{cd}\chi_{ef}^d {\color{blue}R_{bf}^{-1}R_{ae}^{-1}} = \chi_{ab}'^c\underbrace{(R_{ae} R_{bf})(R_{bf}^{-1}R_{ae}^{-1})}_{1}.
 \end{equation}
$$


上述等式表示，三阶张量 $\chi_{ab}^c$ 在$\mathcal{C}_4$对称操作下受到的约束条件是：

$$
\begin{equation}
\chi_{ab}'^c =  R_{cd} \chi_{ef}^dR_{bf}R_{ae}.
 \end{equation}
$$
在得到这个表达式的时候利用了旋转操作$R^{-1}=R^T$的性质。
{:.success}


现在考虑沿$z$方向的磁化在绕$z$轴的转动操作$C_z$下受到约束的结果，此时就有$c=d=z,R_{cd}=1$，因为这个时候无论怎么操作，$z$方向都是不变的，从而得到

$$
\begin{equation}
\chi_{ab}^z = \chi_{ef}^zR_{bf}R_{ae}.
 \end{equation}
$$

$\mathcal{C}_4$操作的表示矩阵为

$$
\begin{equation}
R = 
\begin{pmatrix}
0 & -1 \\
1 & 0
\end{pmatrix}.
 \end{equation}
$$

因此，矩阵元素如下：

$$
\begin{equation}
R_{11} = 0, \, R_{12} = -1,\qquad R_{21} = 1, \, R_{22} = 0
 \end{equation}
$$

将 $R_{ae}$ 和 $R_{bf}$ 代入上面的关系：

$$
\begin{equation}
\chi_{ab}^z = \chi_{ef}^zR_{bf}R_{ae}.
 \end{equation}
$$

我们逐个计算 $\chi_{ab}^z$ 的分量：
### 对角分量
- $\chi_{xx}^z$

$$
\begin{equation}
\chi_{xx}^z = R_{xe} R_{xf} \chi_{ef}^z = R_{11} R_{11} \chi_{xx}^z + R_{11} R_{12} \chi_{xy}^z + R_{12} R_{11} \chi_{yx}^z + R_{12} R_{12} \chi_{yy}^z.
 \end{equation}
$$

代入 $R$ 的元素 $R_{11} = 0$, $R_{12} = -1$：

$$
\begin{equation}
\chi_{xx}^z = (0)(0) \chi_{xx}^z + (0)(-1) \chi_{xy}^z + (-1)(0) \chi_{yx}^z + (-1)(-1) \chi_{yy}^z.
 \end{equation}
$$

简化得到：

$$
\begin{equation}
\chi_{xx}^z = \chi_{yy}^z.
 \end{equation}
$$



- **$\chi_{yy}^z$**

类似地，可以证明 $\chi_{yy}^z = \chi_{xx}^z$。

### 非对角分量

- **$\chi_{xy}^z$ 和 $\chi_{yx}^z$**

对于 $\chi_{xy}^z$ 和 $\chi_{yx}^z$，可以进行类似的计算：

$$
\begin{equation}
\chi_{xy}^z = R_{xe} R_{yf} {\color{blue}\chi_{ef}^z}.
 \end{equation}
$$

代入 ${\color{blue}R_{21} = 1=-R_{12}}$, $R_{22}=R_{11} = 0$，并计算非对角分量(${\color{blue}e=x,f=y}$)，可以得到：

$$
\begin{equation}
\chi_{xy}^z = -\chi_{yx}^z.
 \end{equation}
$$


通过以上计算，$\chi_{ab}^z$ 在 $\mathcal{C}_4$ 旋转对称性下的约束为：

$$
\begin{equation}
\chi_{ab}^z =
\begin{pmatrix}
\chi & \beta \\
-\beta & \chi
\end{pmatrix},
 \end{equation}
$$

其中：
- $\chi = \chi_{xx}^z = \chi_{yy}^z$ 是对角元，表示各向同性分量。
- $\beta = \chi_{xy}^z = -\chi_{yx}^z$ 是非对角元，表示反对称分量。


# 实例讨论
前面讨论的时候是需要知道操作的表示矩阵，但是有时候就只想通过对称性对物理量的约束来分析某一种响应是不是存在，而并不想详细的知道响应系数各个分量之间的依赖关系，这里就整理一下该怎么做。这里就以[Nonlinear Superconducting Magnetoelectric Effect](https://arxiv.org/abs/2404.18616)这篇文章中的结果进行举例，考虑$\mathcal{C}_4\mathcal{T}$对称性对一阶$\delta M_a^{(1)}$以及二阶$\delta M_a^{(2)}$响应的约束。


- $ \mathcal{C}_4 $ 是 **四分之一角旋转操作**，将坐标系旋转 $ 90^\circ $（顺时针或逆时针）。
- 在二维动量空间中，动量分量 $ \mathbf{q} = (q_x, q_y) $ 的旋转由 $ \mathcal{C}_4 $ 作用下的旋转矩阵 $ R $ 描述。

在 $ x $-$ y $ 平面，顺时针 $ 90^\circ $ 旋转的矩阵为：

$$
\begin{equation}
R = 
\begin{bmatrix}
0 & 1 \\
-1 & 0
\end{bmatrix}.
 \end{equation}
$$

对动量 $ \mathbf{q} = (q_x, q_y) $，施加 $ \mathcal{C}_4 $ 操作后的动量变为：

$$
\begin{equation}
\mathbf{q}' = R \mathbf{q} = 
\begin{bmatrix}
0 & 1 \\
-1 & 0
\end{bmatrix}
\begin{bmatrix}
q_x \\
q_y
\end{bmatrix}
=
\begin{bmatrix}
q_y \\
-q_x
\end{bmatrix}.
 \end{equation}
$$

因此，$ \mathcal{C}_4 $ 旋转操作将动量分量 $ (q_x, q_y) $ 变换为：
$
(q_x, q_y) \xrightarrow{\mathcal{C}_4} (q_y, -q_x).
$


- 时间反演操作（记作 $ \mathcal{T} $）会将动量反向，因为时间反演对速度取负号，进而导致动量的反向。

对动量 $ \mathbf{q} = (q_x, q_y) $，时间反演操作的作用为：
$$
\begin{equation}
\mathbf{q} \xrightarrow{\mathcal{T}} -\mathbf{q}.
 \end{equation}
$$

即：

$$
\begin{equation}
(q_x, q_y) \xrightarrow{\mathcal{T}} (-q_x, -q_y).
 \end{equation}
$$


联合操作 $ \mathcal{T} \mathcal{C}_4 $ 是先进行 $ \mathcal{C}_4 $ 旋转操作，再进行时间反演操作。

### 变换关系：
1. **先进行 $ \mathcal{C}_4 $ 旋转**：
  
   $$
   (q_x, q_y) \xrightarrow{\mathcal{C}_4} (q_y, -q_x).
   $$

2. **再进行时间反演 $ \mathcal{T} $**：
   
   $$
   \begin{equation}
   (q_y, -q_x) \xrightarrow{\mathcal{T}} (-q_y, q_x).
   \end{equation}
   $$

因此，联合操作 $ \mathcal{T} \mathcal{C}_4 $ 对动量的变换关系为：

$$
\begin{equation}
(q_x, q_y) \xrightarrow{\mathcal{T} \mathcal{C}_4} (-q_y, q_x).
\end{equation}
$$

因为参考的这篇文献中考虑的是个二维系统，所以Cooper对的动量$\boldsymbol{q}=(q_x,q_y)$，并没有$q_z$分量

- 一阶响应$\delta M_a^{(1)}=\alpha_{ab}q_b,\qquad a={x,y,z},\quad b={x,y}$

首先知道在时间反演对称性$\mathcal{T}$使得二阶响应$\delta M_c^{(2)}$被禁戒，只会存在一阶响应$\delta M_a^{(1)}$，但此时$\mathcal{C}_4$操作对磁化$\delta M_z^{(1)}$是不改变符号的，但是会改变动量

$$
\begin{equation}
\begin{aligned}
&\mathcal{C}_4(\delta M_z^{(1)})\mathcal{C}_4^{-1}=\delta M_z^{(1)}=\alpha_{zb}q_b\quad b={x,y}\\
&b=x:\quad {\color{blue}\mathcal{C}_4(\alpha_{zx}q_x)\mathcal{C}_4^{-1}}=\alpha_{zx}q_x={\color{blue}\alpha_{zx}q_y}\Longrightarrow \alpha_{zx}=0\\
&b=y:\quad {\color{blue}\mathcal{C}_4(\alpha_{zy}q_y)\mathcal{C}_4^{-1}}=\alpha_{zy}q_y={\color{blue}-\alpha_{zy}q_x}\Longrightarrow \alpha_{zy}=0\\
\end{aligned}
\end{equation}
$$

可以发现$\mathcal{C}_z$操作使得一阶响应的$M_z$方向的磁化为零。而$\mathcal{C}_z$操作对于磁化分量的变化为

$$
\begin{equation}
(M_x,M_y,M_z)\xrightarrow{\mathcal{C}_z}(M_y,-M_x,M_z)
\end{equation}
$$

类似的计算对$\delta M_x$和$\delta M_y$分析可以发现$\mathcal{C}_z$对于$x,y$方向的磁化是不会禁戒的。

- 二阶响应$\delta M_c^{(2)}=\chi_{ab}^{c}q_aa_b,\qquad a,b={x,y},\quad c={x,y,z}$


根据最开始的分析，在时间反演$\mathcal{T}$操作下，磁化的二阶响应是被禁戒的，但文献中考虑的系统没有单独的时间反演对称性，而是一个联合的$\mathcal{C}_4\mathcal{T}$，而且此时该对称性使得一阶响应$\delta M_z$被禁戒，但是对于二阶响应可以发现却是可以存在的

$$
\begin{equation}
\begin{aligned}
&\mathcal{C}_4(\delta M_c^{(2)})\mathcal{C}_4^{-1}=\delta M_c^{(2)}=\chi_{ab}^cq_aq_b,\qquad a,b=\{x,y\},\qquad c=\{z\}\\
&a=x,b=y:{\color{blue}\mathcal{C}_4(\delta M_z^{(2)})\mathcal{C}_4^{-1}=-\chi_{xy}^zq_yq_x}=\chi^z_{xy}q_xq_y\Longrightarrow\chi^z_{xy}=0\\
&a=x,b=x:{\color{blue}\mathcal{C}_4(\delta M_z^{(2)})\mathcal{C}_4^{-1}=\chi_{xx}^zq_yq_y}=\chi^z_{xx}q_xq_x\\
&a=y,b=y:{\color{blue}\mathcal{C}_4(\delta M_z^{(2)})\mathcal{C}_4^{-1}=\chi_{yy}^{z}q_xq_x}=\chi^z_{yy}q_yq_y\\
\end{aligned}
\end{equation}
$$

此时可以发现$\chi_{xx},\chi_{yy}$在$\mathcal{C}_4$操作下是非零的，而且前面也已经证明过了这两个系数不仅是非零的，还满足$\chi^z_{xx}=\chi^z_{yy}$。综合前面所有的分析则可得到一个表格
<div align="center">

|          | $\mathcal{P}$ | $\mathcal{T}$ | $\mathcal{PT}$ | $C_4\mathcal{T}$ |
|----------|-----------------|-----------------|------------------|--------------------|
| $M_x$  | $\checkmark(2)$   | $\checkmark(1)$   | $\times$       | $\checkmark(1)$      |
| $M_y$  | $\checkmark(2)$   | $\checkmark(1)$   | $\times$       | $\checkmark(1)$      |
| $M_z$  | $\checkmark(2)$   | $\checkmark(1)$   | $\times$       | $\checkmark(2)$      |

</div>





目前刚接触对称性分析，水平有限，若有错误欢迎指出，我也希望能够学习到。
{:.info}







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

