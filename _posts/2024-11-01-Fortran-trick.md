---
title: Fortran踩坑记录(持续更新中)
tags:  Fortran Code 
layout: article
license: true
toc: true
key: a20241101
pageview: true
# cover: /assets/images/Julia/julia-logo.png
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
最近在用`Fortran`写程序的时候，在意想不到的地方踩了坑，这篇Bolg整理一下。
{:.info}
<!--more-->
# 前言
用`Fortran`写程序主打一个方便，只要把公式翻译成代码就行了，但是在码代码的时候还是遇到了一些意想不到的问题，这里整理一下，方便自己查阅闭坑。

## 费米分布函数

在把费米分布转成代码为
```fortran
function fermi(ek)
    ! 费米分布函数
    integer,parameter::dp = kind(1.0)
    real(dp) fermi,ek,kbt
    kbt = 0.01
    fermi = 1.0/(exp(ek/kbt) + 1)  
    ! fermi = ek 
    return
end
```

但这里有$e^x$指数函数，如果$x$足够大的话，如果是在单精度下面(dp = kind(1.0))，此时分布函数就会溢出。不过就算使用双精度(kind(1.0d0))，还是无法完全避免数据溢出这个问题，所以最好的方式就是给定一个阈值，判断在哪些范围内可以使用费米分布函数，其他范围就是0或者1，取决于该量子态是占据还是非占据
```fortran
function fermi(ek)
    !  万万不能直接用费米分布函数，会存在浮点溢出
    integer,parameter::dp = kind(1.0)
    real(dp) fermi,ek,kbt
    kbt = 0.001
    if(ek/kbt>-40 .and. ek/kbt<40) fermi = 1.0/(exp(ek/kbt) + 1) 
    if(ek/kbt <-40) fermi = 1.0
    if(ek/kbt >40) fermi = 0.0
    return
end 
```
除了上面的截断方式，还可以在能量大于零和小于零的情形下，对费米分布函数做一下变形处理，在数值上避免发散

$$
f(E) = 
\begin{cases}
    \frac{e^{-E/(k_B T)}}{1 + e^{-E/(k_B T)}}, & E > 0 \\
    \frac{1}{1 + e^{E/(k_B T)}}, & E \leq 0
\end{cases}
$$

```fortran
real(dp) function fermi(ek)
    use code_param
    implicit none
    real(dp), intent(in) :: ek
    real(dp) :: beta
    beta = 1.0 / kbt
    if (ek .le. 0.0) then
        fermi = 1.0 / (exp(beta * ek) + 1.0)
    else
        fermi = exp(-beta * ek) / (1.0 + exp(-beta * ek))
    end if
end function fermi
```

如果只是单纯的考虑零温情形，那么就可以使用阶跃函数来代替费米分布
```fortran
function fermi(ek)
    real fermi,ek
    if(ek) < 0 then
      fermi = 1
    else
      fermi = 0
    return
end 
```

有限温度下，对费米分布的导数为

$$
f'(E)=-\frac{1}{4k_BT}sech^2(\frac{E}{2k_BT})
$$

```fortran
real(dp) function fermi_derivative(ek)
    use code_param
    implicit none
    real(dp), intent(in) :: ek
    real(dp) :: beta
    beta = 1.0 / kbt
    fermi_derivative = -0.25 * beta * (1.0 / cosh(0.5 * beta * ek))**2
end function fermi_derivative
```

如果是零温极限，此时分布函数是阶跃函数，那么导数就是$\delta$函数的形式，此时在数值上可以利用高斯分布来代替

$$
f'(E)\approx-\frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{E^2}{2\sigma^2}}
$$

```fortran
real(dp) function fermi_derivative_low_temp(ek)
    ! 极低温下费米分布的导数,这里使用高斯展宽代替
    use code_param
    implicit none
    real(dp), intent(in) :: ek
    fermi_derivative_low_temp = -1.0 / (sqrt(2.0 * pi) * sigma) * exp(-ek**2 / (2.0 * sigma**2))
end function fermi_derivative_low_temp
```
其中的$\sigma$就是展宽，具体这个值取多少需要根据程序计算的曲线平滑程度或者经验进行确定。

## 动态数组赋值

在`Fortran`中经常会用到可变大小的数组，通常会在程序执行过程中来确定数组的大小。但这个时候还是需要对数组进行初始化，但有时候会有错误的初始化，比如

```fortran
real,allocatable::ones(:,:)

allocate(ones(10,10))

do i0 = 1,10
    ones(i0,i0) = 1.0
end do

```
好吧，这是一个错误的示范，因为这里只是对对角元素进行了赋值，但没有对其它位置的元素赋值。所以在后面的程序中如果使用了`ones`这个数组，那么必然就会得到错误的结果，因为其它位置的元素在这里没有赋值，但是`Fortran`对于这种动态类型的数组，不会默认其初始位置都是零，所以那些没有赋值的地方，数据的指向就是奇奇怪怪的地方，计算结果必然是错的。

正确的做法应该是先对确定大小的动态数组的所有元素进行初始化，然后再填入相对应的值
```fortran
real,allocatable::ones(:,:)

allocate(ones(10,10))
ones = 0.0

do i0 = 1,10
    ones(i0,i0) = 1.0
end do

```

- 并行可变大小数组声明以及初始化

在并行计算的时候，声明一个可变大小的数组，在确定其矩阵维度之后要对这个矩阵进行初始化，不然在每个进程中，不计算的部分都具有不同的初始值，最后在收集数据的时候就可能得到错误的数据结果
```fortran
real(dp),allocatable::Veff_vec(:,:)
complex(dp),allocatable::Veff_val(:)
real(dp),allocatable::Veff(:,:)  ! 存储费米面上的相互作用,其维度由费米点的数量决定
real(dp),allocatable::Veff_mpi(:,:)  ! 存储费米面上的相互作用,其维度由费米点的数量决定

if(numk_FS.eq.0)then
    write(*,*)"The Number of Fermi surface points = ",numk_FS
    write(*,*)"----------------------------------------------"
    write(*,*)"                    Error                     "
    write(*,*)"----------------------------------------------"
else
    ! write(*,"(A30,10F20.8)")"Number of Fermi surface points = ",numk_FS

    allocate(Veff_mpi(numk_FS,numk_FS))
    allocate(Veff(numk_FS,numk_FS))
    allocate(Veff_vec(numk_FS,numk_FS))  
    allocate(Veff_val(numk_FS))
    Veff_mpi = 0.0    ! 一定要对这些数据进行初始化
    Veff = 0.0
    Veff_vec = 0.0
    Veff_val = 0.0
end if

nki = floor(indcore * (2.0 * kn + 1)/numcore) - kn
nkf = floor((indcore + 1) * (2.0 * kn + 1)/numcore) - kn - 1
do ikx = nki,nkf
do iky = -kn,kn
!   中间是具体计算过程
enddo
end do

call MPI_Barrier(MPI_COMM_WORLD,ierr)   
call MPI_Reduce(Veff_mpi, Veff, numk_FS**2, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
call MPI_BCAST(Veff, numk_FS**2, MPI_REAL, 0, MPI_COMM_WORLD, ierr) 
```


# 函数声明
一般Fortran默认的都是单精度，可以自行改变成双精度，比如对于零温费米分布函数

```fortran
function fermi(ek)
    ! 零温下的分布函数
    implicit none
    integer,parameter::dp = kind(1.0d0)   ! 双精度
    real(dp), intent(in) :: ek
    real(dp) fermi
    if (ek < 0.0) then
        fermi = 1.0
    else
        fermi = 0.0
    end if
end function fermi
```
在主程序或者其他子过程、函数中调用双精度的费米分布函数时，就需要明确声明精度
```fortran
program main
    implicit none
    integer,parameter::dp = kind(1.0d0)   ! 双精度
    real(dp),external::fermi

    stop
end program main
```
有时候也会习惯在函数声明的时候就直接确定精度
```fortran
real(dp) function fermi(ek)
    ! 零温下的分布函数
    implicit none
    integer,parameter::dp = kind(1.0d0)   ! 双精度
    real(dp), intent(in) :: ek
    if (ek < 0.0) then
        fermi = 1.0
    else
        fermi = 0.0
    end if
end function fermi
```
需要说明的是，对于后面这种方式定义的费米分布，在双精度情形下使用时是**错误**的方式，因为此时在函数定义的时候声明，控制精度的dp是不明确的，为了保险起见还是采用第一种方式来声明函数精度。

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