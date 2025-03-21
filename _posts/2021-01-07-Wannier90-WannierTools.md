---
title: Wannier90产生WannierTools所需的紧束缚模型
tags:  vasp
layout: article
license: true
toc: true
key: a20210107a
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
pageview: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
这里主要是整理一下怎么结合VASP和Wannier90来产生WannierTools所需要的紧束缚模型的数据.
{:.info}
<!--more-->
这里先假定已经熟悉了VASP这个软件,而且服务器上面的VASP编译的时候也已经接口好了Wannier90这个软件,那么在一个计算最后想要进行紧束缚模型的分析,首先要先能从VASP计算得到Wannier90所需要的主要文件,其实也就是一个简单的参数设置就可以
> LWANNIER90 = .TRUE. | .FALSE.

这个值默认是`.FALSE.`,所以如果VASP已经和Wannier90接口好,那么把这个参数变为`.TRUE.`即可,最后就可以得出Wannier90所需要的计算文件
```shell
wannier90.mmn
wannier90.eig
wannier90.amn
wannier90.win
```
得到了上面的四个文件之后,其中`wannier90.win`是主要的控制文件,就和WannierTools中的`wt.in`文件一样,决定你要计算的是什么东西,要怎么计算.

# Wannier90参数设置
因为我还在学习Wannier90,所以具体的内容我也不是特别懂,这里就只是相利用它来计算紧束缚模型哈密顿量,所以这里仅仅学习了哪个控制参数可以用来实现这个功能.我这里使用的是Wannier90的examples中的exalpme02这个例子来通过Wannier90计算它的紧束缚哈密顿量数据.

![png](/assets/images/wannierTools/W90_1.png)

![png](/assets/images/wannierTools/W90_2.png)

有了这些文件之后就可以设置参数进行具体的计算了,先来看看`lead.win`这个文件中的参数,不过我并看不明白每一个参数在做什么,不过为了计算得到紧束缚哈密顿量的数据,需要设置`write_hr`这个参数,它是一个逻辑变量,用来控制是否计算并输出紧束缚哈密顿量的数据(逻辑值)
```shell
! Lead : Tutorial Example 2

 num_wann        =   4
 num_iter        = 20
 write_hr = true  ! 计算并输出紧束缚哈密顿量的数据

! SYSTEM

begin unit_cell_cart
bohr
-4.67775 0.00000 4.67775
 0.00000 4.67775 4.67775
-4.67775 4.67775 0.00000

end unit_cell_cart

begin atoms_frac
Pb 0.00   0.00   0.00
end atoms_frac

begin projections
Pb:sp3
end projections

! KPOINTS

mp_grid : 4 4 4

begin kpoints
0.0000  0.0000   0.0000
0.0000  0.2500   0.0000
0.0000  0.5000   0.0000
0.0000  0.7500   0.0000
0.2500  0.0000   0.0000
0.2500  0.2500   0.0000
0.2500  0.5000   0.0000
0.2500  0.7500   0.0000
0.5000  0.0000   0.0000
0.5000  0.2500   0.0000
0.5000  0.5000   0.0000
0.5000  0.7500   0.0000
0.7500  0.0000   0.0000
0.7500  0.2500   0.0000
0.7500  0.5000   0.0000
0.7500  0.7500   0.0000
0.0000  0.0000   0.2500
0.0000  0.2500   0.2500
0.0000  0.5000   0.2500
0.0000  0.7500   0.2500
0.2500  0.0000   0.2500
0.2500  0.2500   0.2500
0.2500  0.5000   0.2500
0.2500  0.7500   0.2500
0.5000  0.0000   0.2500
0.5000  0.2500   0.2500
0.5000  0.5000   0.2500
0.5000  0.7500   0.2500
0.7500  0.0000   0.2500
0.7500  0.2500   0.2500
0.7500  0.5000   0.2500
0.7500  0.7500   0.2500
0.0000  0.0000   0.5000
0.0000  0.2500   0.5000
0.0000  0.5000   0.5000
0.0000  0.7500   0.5000
0.2500  0.0000   0.5000
0.2500  0.2500   0.5000
0.2500  0.5000   0.5000
0.2500  0.7500   0.5000
0.5000  0.0000   0.5000
0.5000  0.2500   0.5000
0.5000  0.5000   0.5000
0.5000  0.7500   0.5000
0.7500  0.0000   0.5000
0.7500  0.2500   0.5000
0.7500  0.5000   0.5000
0.7500  0.7500   0.5000
0.0000  0.0000   0.7500
0.0000  0.2500   0.7500
0.0000  0.5000   0.7500
0.0000  0.7500   0.7500
0.2500  0.0000   0.7500
0.2500  0.2500   0.7500
0.2500  0.5000   0.7500
0.2500  0.7500   0.7500
0.5000  0.0000   0.7500
0.5000  0.2500   0.7500
0.5000  0.5000   0.7500
0.5000  0.7500   0.7500
0.7500  0.0000   0.7500
0.7500  0.2500   0.7500
0.7500  0.5000   0.7500
0.7500  0.7500   0.7500
end kpoints
```
执行计算
> wannier90.x lead

计算结束之后就可以看到结果输出了

![png](/assets/images/wannierTools/W90_3.png)

得到了`lead_hr.dat`这个文件数据之后,就可以利用WannierTools来计算其相关的性质了.

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