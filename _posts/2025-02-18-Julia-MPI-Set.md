---
title: Mac上配置MPI并行Julia
tags:  Code Julia
layout: article
license: true
toc: true
key: a20250218
pageview: true
cover: /assets/images/Julia/a4.png
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
这里记录一下怎么在本地电脑上配对Julia的MPI库实现并行
{:.info}
<!--more-->
# 安装
Julia下载安装即可，MPI只需要在包安装环境下执行安装

```julia
julia> import Pkg

julia> Pkg.add("MPI")
```

安装完成之后，在用户的主目录下(主目录直接执行*cd*进入即可)会存在一个隐藏的.julia文件夹，这里前面有个英文句号，进入这个文件夹显示如下

![png](/assets/images/Julia/a1.png)

所有安装的包都会整合在**packages**中，进入其中就可以找到安装好的[MPI](https://juliaparallel.org/MPI.jl/v0.20/)

![png](/assets/images/Julia/a2.png)

![png](/assets/images/Julia/a3.png)

进入到MPI文件夹(中间会存在一个奇怪的文件夹，直接进入即可)的bin目录中就可以找到**mpiexejl**这个执行julia进行MPI并行的执行命令了，接下来就是将这个路径在.zshrc中追加到PATH里面就能够使用了

```shell
export PATH="/Users/Yourname/.julia/packages/MPI/TKXAj/bin:$PATH"   # Julia MPIEX
```

添加之后执行`source .zshrc`或者重启终端就能让命令生效了。

# 执行

假如现在有一个使用了[MPI并行的程序](https://yxli8023.github.io/2024/03/24/chi0-mpi.md.html)

```julia
using LinearAlgebra,DelimitedFiles,Printf,MPI,Dates
#-------------------------------------------------------------------------------
function ham(kx::Float64,ky::Float64)
    t0::Float64 = 0.1  # 缩放系数
    t1x::Float64 = -0.483 * t0
    t1z::Float64 = -0.110 * t0
    t2x::Float64 = 0.069 * t0
    t2z::Float64 = -0.017 * t0
    t3xz::Float64 = 0.239 * t0
    t4xz::Float64 = -0.034 * t0
    tvx::Float64 = 0.005 * t0
    tvz::Float64 = -0.635 * t0
    ex::Float64 = 0.776 * t0
    ez::Float64 = 0.409 * t0
    ham = zeros(ComplexF64,4,4)
    ham[1,1] = 2 * t1x * (cos(kx) + cos(ky)) + 4 * t2x*cos(kx)*cos(ky) + ex
    ham[2,2] = 2 * t1z * (cos(kx) + cos(ky)) + 4 * t2z*cos(kx)*cos(ky) + ez
    ham[1,2] = 2 * t3xz * (cos(kx) - cos(ky))
    ham[2,1] = 2 * t3xz * (cos(kx) - cos(ky))

    ham[3,3] = 2 * t1x * (cos(kx) + cos(ky)) + 4 * t2x*cos(kx)*cos(ky) + ex
    ham[4,4] = 2 * t1z * (cos(kx) + cos(ky)) + 4 * t2z*cos(kx)*cos(ky) + ez
    ham[3,4] = 2 * t3xz * (cos(kx) - cos(ky))
    ham[4,3] = 2 * t3xz * (cos(kx) - cos(ky))

    ham[1,3] = tvx
    ham[1,4] = 2 * t4xz * (cos(kx) - cos(ky))
    ham[2,3] = 2 * t4xz * (cos(kx) - cos(ky))
    ham[2,4] = tvz

    ham[3,1] = tvx
    ham[4,1] = 2 * t4xz * (cos(kx) - cos(ky))
    ham[3,2] = 2 * t4xz * (cos(kx) - cos(ky))
    ham[4,2] = tvz
    val,vec = eigen(ham)
    return val,vec
end
#-------------------------------------------------------------------------------
function fsd(x::Float64)
    T::Float64 = 0.001 # Temperature
    return 1/(exp(x/T) + 1)
end
#-------------------------------------------------------------------------------
# 裸的自旋磁化率
function chi0(qx::Float64,qy::Float64,omega::Float64,nk::Int64)
    # nk::Int64 = 100 # 点撒密一点才能找到费米面
    klist = range(-pi,pi,length = nk)
    bearchi0 = zeros(ComplexF64,4,4)
    for kx in klist  
        for ky in klist 
            val,vec = ham(kx,ky)
            valq,vecq = ham(kx + qx,ky + qy)
            for l1 in 1:4,l2 in 1:4
                re1::ComplexF64 = 0
                for m in 1:4,n in 1:4
                    #  计算极化率
                    re1 += (fsd(val[n]) - fsd(valq[m]))/(im * (omega + 0.0001) + val[n] - valq[m]) *
                         vecq[l1,m]' * vecq[l2,m] * vec[l2,n]' * vec[l1,n]
                    # re1 += (fsd(val[n]) - fsd(valq[m]))/(im * (omega + 0.0001) + val[n] - valq[m]) * 
                    #    vecq[l1,m] * vecq[l2,m]' * vec[l2,n] * vec[l1,n]'
                end
                bearchi0[l1,l2] += re1
            end
        end
    end
    return -1/nk^2 * bearchi0
end
#-------------------------------------------------------------------------------
function chi(qx::Float64,qy::Float64,omega::Float64,nk::Int64)
    U0::Float64 = 3.0
    J0::Float64 = 0.4
    a1 = diagm(ones(2))
    a2 = zeros(Float64,2,2)
    I0 = diagm(ones(4))
    a2[1,1] = U0
    a2[2,2] = U0
    a2[1,2] = J0/2
    a2[2,1] = J0/2
    gamma = kron(a1,a2)
    bearchi0 = chi0(qx,qy,omega,nk)
    chitemp = inv(I0 - bearchi0 * gamma) * bearchi0
    #return chitemp
    return imag(sum(chitemp)),real(sum(chitemp))
end
#-------------------------------------------------------------------------------
function main1(nk::Int64)
    # nk::Int64 = 200
    klist = range(0,pi,length = nk)
    qxlist = zeros(Float64,nk^2)
    qylist = zeros(Float64,nk^2)
    chilist = zeros(Float64,nk^2,2)

    #*************************************************
    # Parameter for MPI 
    MPI.Init()
    comm = MPI.COMM_WORLD
    root = 0
    numcore = MPI.Comm_size(comm)  # 申请的核心数量
    indcore = MPI.Comm_rank(comm)  # 核心的id标号
    #*************************************************


    #*************************************************
    # 循环区间分配
    nki = floor(indcore * nk/numcore) + 1
    nkf = floor((indcore + 1) * nk/numcore) 
    if (MPI.Comm_rank(comm) == root)
        println("开始计算极化率: ",Dates.now())
        println("Number of nk : ",nk)
    end
    #*************************************************
    for iqx in nki:nkf
    # for iqx in 1:nk
        for iqy in 1:nk
            i0 = Int((iqx - 1) * nk + iqy)
            iqx = Int(iqx)
            iqy = Int(iqy)
            qxlist[i0] = klist[iqx]
            qylist[i0] = klist[iqy]
            chilist[i0,1],chilist[i0,2] = chi(klist[iqx],klist[iqy],0.0,nk)
        end
    end
    MPI.Barrier(comm)
    qxlist = MPI.Reduce(qxlist,MPI.SUM,root,comm)
    qylist = MPI.Reduce(qylist,MPI.SUM,root,comm)
    chilist = MPI.Reduce(chilist,MPI.SUM,root,comm)

    if (MPI.Comm_rank(comm) == root)
        println("结束计算极化率: ",Dates.now())
        temp1 = (a->(@sprintf "%3.2f" a)).(nk)
        fx1 ="mpi-chi-nk-" * temp1 * ".dat"
        f1 = open(fx1,"w")
        x0 = (a->(@sprintf "%5.3f" a)).(qxlist)
        y0 = (a->(@sprintf "%5.3f" a)).(qylist)
        z0 = (a->(@sprintf "%5.3f" a)).(chilist[:,1])
        z1 = (a->(@sprintf "%5.3f" a)).(chilist[:,2])
        writedlm(f1,[x0 y0 z0 z1],"\t")
        close(f1)
    end
       
end
#-------------------------------------------------------------------------------
# @time hpath()
# @time ferminsurface(1000)
main1(128)
```

执行
```shell
mpiexejl -n 4 julia filename.jl
```
这里就启动了4个核心计算

# Tips
上面的操作对于Linux系统也是同样的操作，只不过需要在.bashrc中执行路径追加操作，Mac则是.zshrc



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

