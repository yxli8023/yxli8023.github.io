---
title: BBH Nested Wilson loop计算(Julia Version)
tags: Julia Topology
layout: article
license: true
toc: true
key: a20220323
pageview: true
# cover: /assets/images/GroupTheory/cube_symmetry.jpg
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
`Python`在执行一些操作的时候有点慢，这里就使用`Julia`写了一下Nested Wilson loop以及Wilson loop的计算程序。
{:.info}
<!--more-->
# Wilson loop
这个计算没有什么好说的,在其他模型里面也已经做过很多次了,比如[BHZ模型Wilson loop计算](https://yxli8023.github.io/2021/06/28/BHZ-Wilson.html)这篇博客,废话不多说直接上代码
```julia
using LinearAlgebra,DelimitedFiles
#--------------------------------------------------------------
function hamset(kx::Float64,ky::Float64)
    hn::Int64 = 4
    gamx::Float64 = 0.5  
    lamx::Float64 = 1.0  
    gamy::Float64 = gamx
    lamy::Float64 = lamx
    xsyb1::Float64 = 0.0000000000000    
    xsyb2::Float64 = 1.0000000000001
    ysyb1::Float64 = 0.0000000000000    
    ysyb2::Float64 = 1.000000000000
    ham = zeros(ComplexF64, hn, hn)
    ham[1, 1] = xsyb1
    ham[2, 2] = ysyb1
    ham[3, 3] = ysyb1
    ham[4, 4] = xsyb1
    ham[1, 3] = (gamx + lamx * exp(im * kx)) * ysyb2
    ham[2, 4] = gamx + lamx * exp(-im * kx)
    ham[1, 4] = gamy + lamy * exp(im * ky)
    ham[2, 3] = (-gamy - lamy * exp(-im * ky)) * xsyb2
    ham[3, 1] = conj(ham[1, 3])
    ham[4, 2] = conj(ham[2, 4])
    ham[4, 1] = conj(ham[1, 4])
    ham[3, 2] = conj(ham[2, 3])
    return ham
end
#--------------------------------------------------------------------------------------
function Wilson_kx(kx::Float64)
    nky::Int64 = 100
    hn::Int64 = 4
    Nocc::Int64 = Int(hn/2)
    wave = zeros(ComplexF64,hn,hn,nky) # 存储哈密顿量对应的波函数
    wan = zeros(ComplexF64,Nocc,Nocc) # 存储Wannier哈密顿量对应的波函数
    F = zeros(ComplexF64,Nocc,Nocc)
    kylist = range(0, 2*pi, length = nky)
    for iy in 1:nky # 固定kx，沿着ky方向计算Wilson loop
        ky = kylist[iy]
        val,vec = eigen(hamset(kx,ky))
        wave[:,:,iy] = vec[:,:] # 存储波函数
    end
    wave[:,:,nky] = wave[:,:,1] # 在边界上波函数首尾相接
    for i1 in 1:Nocc
        F[i1,i1] = 1 # 构建单位矩阵
    end
    for i1 in 1:nky - 1 # index ki lattice
        for i2 in 1:Nocc
            for i3 in 1:Nocc
                wan[i2,i3] = wave[:,i2,i1]' * wave[:,i3,i1 + 1]   # 计算Berry联络
            end
        end
        F = wan * F # 沿着ky方向构造Wannier哈密顿量
    end
    val,vec = eigen(F) # 这里求解得到的本征态是按照本征值的大小进行排列的，和python有一点不同
    # return vec[:,1], vec[:,2]    # 给出两个占据态能带对应的Wannier哈密顿量的本征矢量
    return sort(map(angle,val)/(2*pi))
end
#-------------------------------------------------------------------------------------
function Wilson_ky(ky::Float64)
    nkx::Int64 = 100
    hn::Int64 = 4
    Nocc::Int64 = Int(hn/2)
    wave = zeros(ComplexF64,hn,hn,nkx)
    wan = zeros(ComplexF64,Nocc,Nocc)
    F = zeros(ComplexF64,Nocc,Nocc)
    kxlist = range(0, 2*pi, length = nkx)
    for ix in 1:nkx # 固定ky，沿着kx方向计算Wilson loop
        kx = kxlist[ix]
        val,vec = eigen(hamset(kx,ky))
        wave[:,:,ix] = vec[:,:]
    end
    wave[:,:,nkx] = wave[:,:,1] # 波函数首尾相接
    for i1 in 1:Nocc
        F[i1,i1] = 1
    end
    for i1 in 1:nkx - 1
        for i2 in 1:Nocc
            for i3 in 1:Nocc
                wan[i2,i3] = wave[:,i2,i1]' * wave[:,i3,i1 + 1]   # 计算Berry联络
            end
        end
        F = wan * F # 构造Wannier 哈密顿量
    end
    val,vec = eigen(F)
    # return vec[:,1],vec[:,2]# 给出两个占据态能带对应的Wannier哈密顿量的本征矢量
    return sort(map(angle,val)/(2*pi))
end    
#------------------------------------------------------------
function Wilson(num)
    nkx = 100
    x1list = zeros(Float64,nkx,2)
    x2list = zeros(Float64,nkx,2)
    kxlist = range(-pi,pi,length = nkx)
    i0 = 1
    for kx in kxlist
        append!(kxlist,kx)
        x1 = Wilson_kx(kx)
        x2 = Wilson_ky(kx)
        x1list[i0,:] = x1
        x2list[i0,:] = x2
        i0 += 1
    end
    fx1 = "wilson-kx-" * string(num) * ".dat"
    f1 = open(fx1,"w")
    writedlm(f1,[kxlist x1list])
    close(f1)

    fx1 = "wilson-ky-" * string(num) * ".dat"
    f1 = open(fx1,"w")
    writedlm(f1,[kxlist x2list])
    close(f1)
end
#-----------------------------------------------------------------
@time Wilson(1)
```

# Nested Wilson Loop
关于Nested Wilson loop的计算可以参考[Electric multipole moments, topological multipole moment pumping, and chiral hinge states in crystalline insulators
](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.245115)这篇文章，这里主要说一下自己在学习计算的时候踩过的坑。

与计算Wilson loop相同，这里最主要的仍然是找到一个Wannier band basis，也就是文章的中的公式

$$\rvert w_{x,\mathbf{k}}^j\rangle=\sum_{n=1}^\text{Nocc}\rvert u^n_\mathbf{k}\rangle[v_{x,\mathbf{k}}^j]^n$$

其实在做计算的时候，最让人困扰的不过是公式中的一大堆符号对应的到底是什么，这里就来讲这个公式拆解开，一步一步的讲清楚里面的含义。这里假设你已经知道为什么要计算Nested Wilson loop，
我在这里就简单的阐述一下。首先要是体系的Wilson loop计算对应的Wannier哈密顿量的能带是有能隙的，也就是说你的体系是4带模型，那么当占据态是2条能带的时候，每个占据态能带会对应着一个Wannier center，
比如BHZ模型的两条占据态能带对应的Wannier band就是相互交叉的，而且因为Wilson loop与边界态之间的拓扑等价性，TI是有边界态的，所以其对应的Wilson loop在形状上就与边界态类似。而对于高阶拓扑相，
首先就是要使得边界态打开能隙，那么相对应的就是其Wilson loop计算得到的Wannier center随着某个动量参数的演化是不会相互交叉的，这一点在上面BBH模型中已经计算过了，所以此时就可以对某一个单独的Wannier band
计算它的Nested Wilson loop，所以首先第一步就是必须要明白什么样的情况下，是需要计算体系的Nested Wilson loop。

这里的$\sum_{n=1}^\text{Nocc}$不用讲太多，是需要对占据态进行求和，但是这个$n$其实表示的只是说哈密顿量的占据态，也就是说对于$\rvert u^n_\mathbf{k}\rangle$而言，这是哈密顿量的占据态波函数，$n$表示占据态其实是对$\rvert u^n_\mathbf{k}\rangle$
而言的，虽然$[v_{x,\mathbf{k}}^j]^n$中同样存在这个$n$，但是在这个地方$n$代表的不是占据态，在这里$j$代表的才是你选择的是哪一个Wannier band来计算其对应的Nested Wilson loop，也就是这里$j$代表的是你选择的是那个占据的Wannier band，而$n$在这里
表示的是一个Wannier哈密顿量本征矢量的第$n$个分量。假如$H_\text{Wann}$是Wannier哈密顿量，其满足


$$H_\text{Wann}\rvert v_\mathbf{k}\rangle=E\rvert v_\mathbf{k}\rangle$$

那么这里的$[v_{x,\mathbf{k}}^j]^n$表示的就是这个本征矢量的第$n$个分量，$j$则表示你选择是哪个本征值对应的本征矢量，也就是选择了哪一个Wannier band。这里的$x$则表示你在做Wilson loop的时候，是沿着哪个方向进行的，即就是讲上面公式中的$H_\text{Wann}$替换成你
构建的那个Wilson loop的哈密顿量就可以。

至于$\rvert u^n_\mathbf{k}\rangle$就很简单了，它表示的就是你的哈密顿量的本征态，当然了在计算的时候，还是要选择正确的占据态才可以。下面直接上代码，在其中同样做了注释
```julia
using LinearAlgebra,DelimitedFiles
#--------------------------------------------------------------
function hamset(kx::Float64,ky::Float64)
    hn::Int64 = 4
    gamx::Float64 = 0.5  
    lamx::Float64 = 1.0  
    gamy::Float64 = gamx
    lamy::Float64 = lamx
    xsyb1::Float64 = 0.000000000000    
    xsyb2::Float64 = 1.1
    ysyb1::Float64 = 0.000000000000    
    ysyb2::Float64 = 1.000000000000
    ham = zeros(ComplexF64, hn, hn)
    ham[1, 1] = xsyb1
    ham[2, 2] = ysyb1
    ham[3, 3] = ysyb1
    ham[4, 4] = xsyb1
    ham[1, 3] = (gamx + lamx * exp(im * kx)) * ysyb2
    ham[2, 4] = gamx + lamx * exp(-im * kx)
    ham[1, 4] = gamy + lamy * exp(im * ky)
    ham[2, 3] = (-gamy - lamy * exp(-im * ky)) * xsyb2
    ham[3, 1] = conj(ham[1, 3])
    ham[4, 2] = conj(ham[2, 4])
    ham[4, 1] = conj(ham[1, 4])
    ham[3, 2] = conj(ham[2, 3])
    return ham
end
#--------------------------------------------------------------------------------------
function Wilson_kx(kx::Float64,nky::Int64)
    # nky::Int64 = 100
    hn::Int64 = 4
    Nocc::Int64 = Int(hn/2)
    wave = zeros(ComplexF64,hn,hn,nky) # 存储哈密顿量对应的波函数
    wan = zeros(ComplexF64,Nocc,Nocc) # 存储Wannier哈密顿量对应的波函数
    F = zeros(ComplexF64,Nocc,Nocc)
    kylist = range(0, 2*pi, length = nky)
    for iy in 1:nky # 固定kx，沿着ky方向计算Wilson loop
        ky = kylist[iy]
        val,vec = eigen(hamset(kx,ky))
        wave[:,:,iy] = vec[:,:] # 存储波函数
    end
    wave[:,:,nky] = wave[:,:,1] # 在边界上波函数首尾相接
    for i1 in 1:Nocc
        F[i1,i1] = 1 # 构建单位矩阵
    end
    for i1 in 1:nky - 1 # index ki lattice
        for i2 in 1:Nocc
            for i3 in 1:Nocc
                wan[i2,i3] = wave[:,i2,i1]' * wave[:,i3,i1 + 1]   # 计算Berry联络
            end
        end
        F = wan * F # 沿着ky方向构造Wannier哈密顿量
    end
    val,vec = eigen(F) # 这里求解得到的本征态是按照本征值的大小进行排列的，和python有一点不同
    vec1 = vec[:,sortperm(map(angle,val))[1]]
    vec2 = vec[:,sortperm(map(angle,val))[2]]
    return vec1,vec2# 给出两个占据态能带对应的Wannier哈密顿量的本征矢量
end
#-------------------------------------------------------------------------------------
function Wilson_ky(ky::Float64,nkx::Int64)
    # nkx::Int64 = 100
    hn::Int64 = 4
    Nocc::Int64 = Int(hn/2)
    wave = zeros(ComplexF64,hn,hn,nkx)
    wan = zeros(ComplexF64,Nocc,Nocc)
    F = zeros(ComplexF64,Nocc,Nocc)
    kxlist = range(0, 2*pi, length = nkx)
    for ix in 1:nkx # 固定ky，沿着kx方向计算Wilson loop
        kx = kxlist[ix]
        val,vec = eigen(hamset(kx,ky))
        wave[:,:,ix] = vec[:,:]
    end
    wave[:,:,nkx] = wave[:,:,1] # 波函数首尾相接
    for i1 in 1:Nocc
        F[i1,i1] = 1
    end
    for i1 in 1:nkx - 1
        for i2 in 1:Nocc
            for i3 in 1:Nocc
                wan[i2,i3] = wave[:,i2,i1]' * wave[:,i3,i1 + 1]   # 计算Berry联络
            end
        end
        F = wan * F # 构造Wannier 哈密顿量
    end
    val,vec = eigen(F)
    vec1 = vec[:,sortperm(map(angle,val))[1]]
    vec2 = vec[:,sortperm(map(angle,val))[2]]
    return vec1,vec2# 给出两个占据态能带对应的Wannier哈密顿量的本征矢量
    # return sort(map(angle,val)/(2*pi))
end
#------------------------------------------------------------------------------------
function Nseted_Wilson_loop_kx(nkx::Int64)
    # nkx::Int64 = 100
    nky::Int64 = nkx
    hn::Int64 = 4
    Nocc::Int64 = Int(hn/2)
    kxlist = range(-pi, pi, length = nkx)
    kylist = range(-pi, pi, length = nky)
    wave = zeros(ComplexF64,hn,hn,nky)
    pmulist = []
    for ix in 1:nkx
        kx = kxlist[ix]
        for iy in 1:nky
            ky = kylist[iy]
            val,vec = eigen(hamset(kx,ky)) # 计算哈密顿量对应的本征矢量
            wave[:,:,iy] = vec[:,:]
        end
        wave[:,:,nky] = wave[:,:,1] # 波函数首尾相接

        wmu = zeros(ComplexF64,hn,nky)  # 用来构建新的Wannier basis
        for iy in 1:nky
            ky = kylist[iy]
            wann_v1, wann_v2 = Wilson_ky(ky,nkx) # 在固定ky的情况下，计算沿着kx方向的Wilson loop并得到对应的本征矢量
            wmu[:,iy] = wave[:,1,iy] * wann_v1[1] + wave[:,2,iy] * wann_v1[2] # 构建新的Wannier basis
        end
        wmu[:,nky] = wmu[:,1] # 首尾相接

        # 在新的形式下构建Wilson loop
        wan = 1
        for iy in 1:nky - 1
             F0 = wmu[:, iy]' * wmu[:, iy + 1] # 在新的Wannier basis下面构建Wilson loop，也就是计算Nested Wilson loop
             wan = F0 * wan
        end
        pmu = log(wan)/(2 * im * pi)
        if real(pmu) < 0
            pmu += 1
        end
        append!(pmulist,real(pmu))
    end
    return kxlist,pmulist
end
#-----------------------------------------------------------------------
function Nseted_Wilson_loop_ky(nkx::Int64)
    # nkx = 100
    nky = nkx
    hn = 4
    Nocc = Int(hn/2)
    kxlist = range(-pi,pi,length = nkx)
    kylist = range(-pi,pi,length = nky)
    wave = zeros(ComplexF64,hn,hn,nky)
    pmulist = []
    for iy in 1:nky
        ky = kylist[iy]
        for ix in 1:nkx
            kx = kxlist[ix]
            val,vec = eigen(hamset(kx,ky)) # 计算哈密顿量对应的本征矢量
            wave[:,:,ix] = vec[:,:]
        end
        wave[:,:,nkx] = wave[:,:,1]

        wmu = zeros(ComplexF64,hn,nkx)
        for ix in 1:nkx
            kx = kxlist[ix]
            wann_v1, wann_v2 = Wilson_kx(kx,nkx) # 在固定ky的情况下，计算沿着kx方向的Wilson loop并得到对应的本征矢量
            wmu[:,ix] = wave[:,1,ix] * wann_v1[1] + wave[:,2,ix] * wann_v1[2] # 构建新的Wannier basis
        end
        wmu[:,nkx] = wmu[:,1] # 首尾相接
        # 在新的形式下构建Wilson loop
        wan = 1
        for ix in 1:nkx - 1
             F0 = wmu[:,ix]' * wmu[:,ix + 1] # 在新的Wannier basis下面构建Wilson loop，也就是计算Nested Wilson loop
             wan = F0 * wan
        end
        pmu = log(wan)/(2*im*pi)
        if real(pmu) < 0
            pmu += 1
        end
        append!(pmulist,real(pmu))
    end
    return kylist,pmulist
end
#----------------------------------------------------------------------------------------------
function  Nested(num)
    nkx = 100
    x1,x2 = Nseted_Wilson_loop_ky(nkx)
    fx1 = "ky-" * string(num) * ".dat"
    f1 = open(fx1,"w")
    writedlm(f1,[x1 x2],"\t")
    close(f1)

    x1,x2 = Nseted_Wilson_loop_kx(nkx)
    fx1 = "kx-" * string(num) * ".dat"
    f1 = open(fx1,"w")
    writedlm(f1,[x1 x2],"\t")
    close(f1)
end
#-----------------------------------------------------------------
@time Nested(1)
```
这里说一下自己在用`julia`编写的时候遇到的一个小问题，因为这个模型想要有well defined的Wilson loop，就必须破坏一下对称性让极化有一个确定的符号选择(具体我在说什么可以去看参考中的这篇文章)。所以就需要在哈密顿量中参数设置的时候**xsyb2::Float64 = 1.1**让他不是1就可以了，这样算出来的结果才是具有确定符号的，否则就会是乱跳的量。
{:.warning}

# 参考
- 1.[Electric multipole moments, topological multipole moment pumping, and chiral hinge states in crystalline insulators
](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.245115)

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