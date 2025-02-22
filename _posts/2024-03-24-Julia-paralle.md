---
title: Julia中MPI,@distributed,@threads三种并行方法的比较
tags:  Julia MPI 
layout: article
license: true
toc: true
key: a20240324
pageview: true
cover: /assets/images/Julia/julia-logo.png
header:
  theme: dark
  background: 'linear-gradient(to right, rgb(101, 78, 163), rgb(234, 175, 200))'
article_header:
  type: overlay
  theme: dark
  background_color: false
  background_image: 
    gradient: 'linear-gradient(to right, #d3cce3, #e9e4f0)'
    image: false
aside:
    toc: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
最近折腾了一下Julia的并行操作，目前我了解到的有三种：MPI并行、@distributed、@threads，这里就用极化率的计算来对比一下三者的速度。
{:.info}
<!--more-->
这里先单独列举一下三种并行方式在语法上的不同
# MPI
```julia
using LinearAlgebra,DelimitedFiles,Printf,MPI,Dates
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
```
# @distributed
```julia
@everywhere function main1(nk::Int64)
    # nk::Int64 = 200
    klist = range(-pi,pi,length = nk)
    qxlist = SharedArray(zeros(Float64,nk^2))
    qylist = SharedArray(zeros(Float64,nk^2))
    chilist = SharedArray(zeros(Float64,nk^2,2))
    @sync @distributed for iqx in 1:nk
        for iqy in 1:nk
            i0 = (iqx - 1) * nk + iqy
            qxlist[i0] = klist[iqx]
            qylist[i0] = klist[iqy]
            chilist[i0,1],chilist[i0,2] = chi(klist[iqx],klist[iqy],0.0,nk)
        end
    end
    temp1 = (a->(@sprintf "%3.2f" a)).(nk)
    fx1 ="distributed-chi-nk-" * temp1 * ".dat"
    f1 = open(fx1,"w")
    x0 = (a->(@sprintf "%5.3f" a)).(qxlist)
    y0 = (a->(@sprintf "%5.3f" a)).(qylist)
    z0 = (a->(@sprintf "%5.3f" a)).(chilist[:,1])
    z1 = (a->(@sprintf "%5.3f" a)).(chilist[:,2])
    writedlm(f1,[x0 y0 z0 z1],"\t")
    close(f1)
end
```
在使用@distributed这个宏的时候，记得要在所有的函数前面加上@everywhere这个宏，让所有的核心都知道要用到的函数。

# @threads
```julia
using SharedArrays, LinearAlgebra,Distributed,DelimitedFiles,Printf,.Threads
#-----------------------------------------------------------------------------
function main1(nk::Int64)
    # nk::Int64 = 200
    klist = range(-pi,pi,length = nk)
    qxlist = SharedArray(zeros(Float64,nk^2))
    qylist = SharedArray(zeros(Float64,nk^2))
    chilist = SharedArray(zeros(Float64,nk^2,2))
    @threads for iqx in 1:nk
        for iqy in 1:nk
            i0 = (iqx - 1) * nk + iqy
            qxlist[i0] = klist[iqx]
            qylist[i0] = klist[iqy]
            chilist[i0,1],chilist[i0,2] = chi(klist[iqx],klist[iqy],0.0,nk)
        end
    end
    temp1 = (a->(@sprintf "%3.2f" a)).(nk)
    fx1 ="threads-schi-nk-" * temp1 * ".dat"
    f1 = open(fx1,"w")
    x0 = (a->(@sprintf "%5.3f" a)).(qxlist)
    y0 = (a->(@sprintf "%5.3f" a)).(qylist)
    z0 = (a->(@sprintf "%5.3f" a)).(chilist[:,1])
    z1 = (a->(@sprintf "%5.3f" a)).(chilist[:,2])
    writedlm(f1,[x0 y0 z0 z1],"\t")
    close(f1)
end
```
在使用@threads的时候，和@distributed方法不一样，此时不需要在每个函数前面加上@everywhere，但是二者都使用了SharedArray这个函数，可以让每个过程同时对一个数组进行幅值操作且不发生互斥。这里还要强调一下@threads这个宏需要导入`.Threads`这个包，前面这里是有一个句号的，或者可以直接使用`Threads.@threads`,这样就不同在`using`的时候导入这个包了。

# 完整代码
需要使用得到的其他函数如下，只要把上面的不同并行方法套用进来就行，不过要注意前面提到的加@everywhere等问题。
```julia
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
    bearchi0 = SharedArray(zeros(ComplexF64,4,4))
    #@sync @distributed for kx in klist
    for kx in klist
        for ky in klist 
            val,vec = ham(kx,ky)
            valq,vecq = ham(kx + qx,ky + qy)
            for l1 in 1:4,l2 in 1:4
                re1::ComplexF64 = 0
                for m in 1:4,n in 1:4
                    re1 += (fsd(val[n]) - fsd(valq[m]))/(im * (omega + 0.0001) + val[n] - valq[m]) * vecq[l1,m]' * vecq[l2,m] * vec[l2,n]' * vec[l1,n]
                    # re1 += (fsd(val[n]) - fsd(valq[m]))/(im * (omega + 0.0001) + val[n] - valq[m]) * vecq[l1,m] * vecq[l2,m]' * vec[l2,n] * vec[l1,n]'
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
```
# 结果对比
计算的时候选择`nk=40`，使用的核数也是40，运行时间如下
- MPI
```shell
======== Job starts at 2024-03-24 14:47:33 on n33 ======== 
开始计算极化率: 2024-03-24T14:47:40.915
Number of nk : 40
结束计算极化率: 2024-03-24T14:47:52.785
======== Job ends at 2024-03-24 14:47:53 on n33 ========
```
- @distributed
```shell
======== Job starts at 2024-03-24 15:35:28 on n40 ======== 
 65.540391 seconds (1.57 M allocations: 103.933 MiB, 0.07% gc time, 2.65% compilation time)
======== Job ends at 2024-03-24 15:36:46 on n40 ========
```
- @threads
```shell
======== Job starts at 2024-03-24 15:34:24 on n12 ======== 
108.048878 seconds (11.88 G allocations: 328.220 GiB, 50.20% gc time, 134.02% compilation time)
======== Job ends at 2024-03-24 15:36:16 on n12 ========
```
从运行结果上来看使用MPI的效率是最高的，但这里使用`@threads`的运行时间居然变得很长，我怀疑应该是我使用的方法有问题，但是其运行效率应该也不会比MPI快。

提交任务的脚本如下
```shell
#!/bin/bash
#SBATCH --job-name=chi0val
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH -t 3-00:00:00

# load the environment
module purge
module load compiler/intel/2021.3.0
module load mpi/intelmpi/2021.3.0
module load apps/vasp/6.4.0-i21
module load apps/wannier90/3.1-i21
#/soft/vasp6_d4/bin

# where is your binary file
source /soft/anaconda/anaconda3/etc/profile.d/conda.sh
julia=/soft/julia-1.9.4/bin/julia
python=/soft/anaconda/anaconda3/bin/python3.11

NUM_MPI=$SLURM_NTASKS


echo "======== Job starts at `date +'%Y-%m-%d %T'` on `hostname` ======== "

#mpirun -np ${NUM_MPI} julia mpi-rpa.jl
#mpirun -np ${NUM_MPI} julia threads-rpa.jl
#julia -t ${NUM_MPI} threads-rpa.jl
julia -p ${NUM_MPI}  distributed-rpa.jl
# python 01-graphene_prim_model.py  

echo "======== Job ends at `date +'%Y-%m-%d %T'` on `hostname` ======== "
```

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