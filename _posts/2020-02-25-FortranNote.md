---
title: Fortran使用中笔记
tags: Code Fortran
layout: article
license: true
key: a20200225
toc: true
pageview: true
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
将看到的Fortran中常用的一些语法记录起来，方便查阅。
<!--more-->

# Fortran中常用的库函数总结

## 1. Cheevd

> 厄密矩阵本征值与本征矢的求解

## 2. getrf,getri 

> getrf 对一个矩阵进行LU分解

> getri 计算由LU分解后矩阵的逆 

# Fortran中一些混乱概念的收集

## 1.Fortran中格式化输出总结

> write(*,'(I3)')num 表示以3个字符宽度来在屏幕上显示num，若num是一个长度小于等于3的内容，均可以正确显示，当它的长度大于3时，则无法正确显示其结果，屏幕上会输出**号，表示无法正确的显示

### 整数输出(I)

整数的输出用字母I来控制。通用格式为 rIw，其中r代表重复的次数，I代表整数(Integer)，w控制表示的长度。

> 例如 3I10:3个整数，每个整数的显示长度是10个宽度单位，他们会在同一行中显示(若是在写入文件的时候)。

### 浮点数输出(F)

 浮点数的控制方法和整数相同，唯一多的控制内容是小数点后面的显示位数  **rFw.d**

在这里d表示小数点后面用多少位来显示，它的长度是包含在w中的

> 例如:write(*,'(1F6.3)')3.5，每个浮点数控制长度为6，小数点后面显示3位，而对于小数点前面的整数部分，则由空格来补足其长度  

### 科学计数(E)

E是用来显示科学计数的方式 **rEw.d**，一般情况下科学计数法下宽度描述必须满足*w>=d + 7*

该控制并不是我们常用的科学计数法显示，他会将数据控制显示为**aEb**的形式，此时a被归一化为0.1-->1.0 之间的某个数，b则代表的是该数据幂指数的次数

### 科学计数(ES)

该控制方法才是我们通常所用的科学计数法显示，它与上面的输出控制性质完全相同除了此时a是一个1.0–>10之间的一个数，它同样需要满足**w>=d+7**的限制条件，否则会显示为一些*号。

### 逻辑输出(L)

控制输出都具有相同的控制方式**rLw**，对于逻辑值，它只有T和F的显示，分别代表真和假

### 字符输出(A)

字符控制由两种可以选择的控制方式:1.**rA**,2.**rAw**

其中第一种是该字符串有多长，则它的控制长度就是多长。第二种则是要求了该字符串的显示长度是w，当字符串长度小于w则以右对齐的方式通过空格来占位，当字符串长度大于w的时候，则只能显示w长度的字符串，多余的部分则不会显示。

### 空格占位符(X)

nX表示在输出行上面增加n个空格，以此来达到控制输出样式的目的。

### 格式描述符的重复执行

从前面可以看到，一个格式控制是可以通过在前面增加数字，来使其重复执行，其实对于两个不同的控制符，通过可以将其组合在一起，是他们也可以重复的执行。

> 2(I6,2F6.3)这样的组合形式也是可以的，表示有6个数，他们两两分组显示为(I6,2F6.3)

### 换行输出符(/)

斜线(/)表示要换行输出

> write(*,'(F6.2/,I6)')-12.3,5虽然表达式上-12.3和5是在同一行中，但是由于/的作用，5将会在下一行的靠头显示

## 2.Fortran中的文件读写操作

### open中常用的参数选项

> unit = number 同一个整数来将打开的文件和这个编号相关联
>
> file = character 该参数用来指定打开文件的文件名(参数必须是个字符串形式)
>
> status = character 用来确定文件的状态，选项有：old，new，replace，scratch，unknown；当指定文件状态为new时，表示该文件之前不存在，程序会创建一个新的文件。使用replace状态时，则会新建一个文件，将原有的文件信息进行覆盖。scratch状态是在程序运行过程中创建一个临时的文件，一旦程序运行结束或者终止，该文件就会消失，不会作为结果文件存在。
>
> action = read/write/readwrite 确定打开文件要执行的操作
>
> iostat = integer number 这个参数将会返回一个整数值，通过可以确定文件是否成功的打开。若返回值为0则代表文件打开成功，返回其他的值代表文件打开失败
>
> iomsg = character 这个参数会返回文件打开状态的信息。若文件成功的打开，则character变量信息不变，若文件打开失败，则会在character中返回文件打开失败的信息

### read 和 write中的iostat与iomsg

read 和 write中，iostat会接受一个返回值，当它为0是表示文件可以进行正常的读写操作，此时iomsg也会返回一个字符串，它与原来该字符串的值一致。若iostat返回一个非0值的时候，表示文件不能进行正常的读写操作，这时候在iomsg中则会返回产生这中错误读写的原因。**这个附加功能在程序中有文件读写的操作时，是个很有用的返回信息**。有这样的信息返回的话，如果出现文件的读写操作，程序不会异常中断，而是将相应的错误信息返回的所定义的变量中，后面的程序还是可以正确的运行的。

## 3.Fortran中的数组

### Fortran中常用数组赋值方法

fortran中，对于小数组有两种赋值方式：

> integer,dimension(5)::a
>
> a = (/ 1,2,3,4,5 /)
>
> a = [1,2,3,4,5]

还可以是用隐式循环来进行数组的赋值

> a = (/ (i, i = 1,5) /)

### Fortran中数组索引的修改

我们还可以在数组声明的时候,直接对数组可用的索引进行声明

> integer,dimension(5)::a,可用索引为1,2,3,4,5
>
> integer,dimension(-2:2)::a,可用索引为-2,-1,0,1,2
>
> integer,dimension(5:9)::a,可用索引为5,6,7,8,9

### 通过常数声明方便控制数组的大小

在编写程序的过程中,通常会遇到要明确数组的大小,但是在程序的调试过程中一般可以选取小的数组来进行调整来验证程序的结果正确与否,而在正式运行程序的时候,通常这个数组的都是很大的,这就要求我们可以方便的通过一个常数声明来便捷的修改数组的大小:

> integer,parameter::nsize = 100
>
> real,dimension(nsiza)::a1
>
> real,dimension(nsize)::a2

通过这样的方式,我们可以在不同的时候,仅仅通过nsize的修改就可以完全的改变程序中数组的大小,从而避免了逐一修改的麻烦.

### 数组元素的引用

数组元素除了可以通过单个的整数索引来读取数组中的值之外,还可以通过另外一个数组作为它的索引,来获取这些索引所对应的数组中的值

> integer,dimension(3)::vec = (/ 1,2,1 /)
>
> integer,dimension(5)::a = [1,2,3,4,5]
>
> a(vec) 这样的数组索引是可以被接受的,他表示分别取出a(1),a(2),a(1)

### 数组的格式化输入与输出

在Fortran中的格式化输出中,r表示该输出样式的重复次数,这一控制在输出数组是显得非常方便.

> integer,dimension(5)::a = [1,2,3,4,5]
>
> write(*,'(5I3)')a,这样就可以将整个数组完全打印到屏幕上

同样也可以使用do的隐式循环来输出数组

> integer,dimension(5)::a = [1,2,3,4,5]
>
> write(*,'(5I3)')(a(i),i = 1,5)

### do隐式循环的嵌套使用

上面隐式循环的使用,是针对一维数组进行的,然而当数组是二维或者更高维的时候,我们则需要使用嵌套的do隐式循环来对数组进行输出

> integer,dimension(3,3)::a
>
> write(*,'(9I4)')((a(i,i),i = 1,3),j = 1,3)

## 4.Fortran 中的过程

### 子过程中的Intent属性

在构建子过程时，可以对其参数列表进行属性设置，明确哪些变量时输入变量，哪些变量时输出变量，这时候需要通过Intent进行变量的属性说明

> intent(in)    形参只用于向子过程传递输入数据
>
> intent(out)	形参仅用于将结果返回给调用程序
>
> intent(inout)	形参同时具有 in 和out的两种属性

在子过程中一旦对每个形参的属性都进行了声明之后，如果对这些参数村咋不符合其性质的操作，那么程序就会报错。

在子过程中，建议给所有的形参据设定其属性，这是因为Fortran中的变量传递方式是**传址调用**，如果不进行特定的属性声明的话，有可能会对传入的参数进行了修改，从而使程序的运行得到错误的结果，而对形参设定了确定的读写性质后，这种不正确的修改则会使程序报错，从而省去了一些调试的麻烦。

### Fortran中用模块来进行共享数据

如果某个子程序或者函数中通过Use来使用某个模块，那么这个程序单元就可以使用模块中定义好的数据变量。

> module share
>
> implicit none
>
> save
>
> integer::a
>
> end module share

在定义模块的时候,加入save语句后,可以保证模块中的变量只是进行共享,而不能对其进行修改的操作,一旦对这些通过模块共享的变量进行修改,程序就会报错.

### 模块过程

> module share
>
> contains
>
> ​	subroutine sub1(i)
>
> ​	implicit none
>
> ​	integer(kind=8),intent(in)::i
>
> ​	………………………..
>
> ​	end subroutine sub1
>
> end module share

如果将一个子程序写道模块内，当我们在程序中使用这个模块并且调用这个模块中的子程序时,编译器可以检查我们在调用这个子过程的时候,给它输入的实际参数类型和它在定义时的类型是否匹配,如果他们的类型不匹配则编译器将会报错.例如上面的share模块,它包含了一个接受整型变量的子过程sub1,如果我们在使用sub1的时候给它输入一个浮点类型的参数:sub1(0.3)那么程序在编译执行的过程中就会报错.

相反,直接定义咋程序中的子过程,即使我们有对子过程形参类型的声明,当我们在程序中电泳它时,我们输入错误的参数类型,程序仍然可以正常的运行,只是运行的结果是不正确的.  而且重要的是,即使我们此时调用子过程的时候,输入的参数数量比定义子过程的时候还要多,那么程序也不会报错,还是可以正常的运行,只不过此时的结果是错误的.



## 5.Fortran中的函数

因为Fortran中参数的传递类型是**传址调用**,所以在函数的调用时一定要保证形参不会出现在左侧,否则就会修改了传入变量的值,而为了避免这种情况,需要在函数定义的时候,对形参使用intent(in)的属性设置,来避免无意中对输入参数值的改变.

## 6.过程作为参数传递给其它过程

通常在程序中，传递给子过程或者函数的参数，都是传递了变量的内存地址。在程序中，函数和子过程，你同样可以将它们看作是“变量”，只不过这种变量透着高级的结构而已，所以从这个角度来看，同样可以将函数和子过程当作参数来进行传递。

### 自定义函数作为参数传递

> subroutine test(func,a,b,result)
>
> implicit none
>
> **real,eternal::func**  !在这里通过声明确定func是个函数,且该函数的结果是个实数
>
> real,intent(in) a,b
>
> real,intent(out) result
>
> result = func(a,b)
>
> end subroutine test

如上所示,我们定义了一个子过程,它可以接受函数作为参数,所以在形参类型声明的时候,用关键字external来确定func是个外部定义的函数.通过这样的定义,在使用该子过程的时候,我们可以将任意一个两个参数输入的函数作为一个"变量"来传递给该子过程,从而得到结果.

### 子程序作为参数传递

> subroutine test(sub1,a,b,result)
>
> implicit none
>
> **eternal::sub1**  !在这里通过声明确定func是个函数,且该函数的结果是个实数
>
> real,intent(in) a,b
>
> real,intent(out) result
>
> result = func(a,b)
>
> end subroutine test

与函数相同,子程序也可以做为参数进行传递,不过与函数不同的一点是,在声明子过程的时候,不需要去定义子过程的返回类型.函数是有确定的返回类型的,而子过程则没有返回类型,它的结果一般都是通过参数,来进行传递的.

## 7.数组的高级特性

> integer,dimension(12,12)::a
>
> integer,dimension(-12:12,-12:12)::b

二维数组的声明和一维数组相似,同样可以直接通过整数来确定数组的大小.同样的我们也可以自己定义数组的索引,从而构建数组.

数组是内存中的一串连续的单元,在Fortran中数组数据的存储时按列存储的.

reshape函数,可以将数组行列大小进行重新调整.reshape函数中,第一个参数是需要改变的数组,第二个参数是希望矩阵改编成什么形状,函数的返回值则是结果.

> program ex01
> implicit none
> integer,dimension(10)::a1 = [1,2,3,4,5,6,7,8,9,10]
> integer,dimension(2,5)::a
> integer i,j
> a = reshape(a1,(/2,5/))
> write(*,*)a1
> write(*,*)a
> do i = 1,2
> 	do j = 1,5
> 		write(*,*)a(i,j)
> 		end do
> 	end do
> stop
> end

从上面这个例子中可以看到,将一个一维的数组通过reshape函数改编成一个2维的矩阵.

同样的,在这里需要注意一个问题,如果要用确定的一连串数值对一个二维的数组进行初始化,那么就需要用到reshape函数来对这个数据进行维数改变,从而让它成为二维数组的初始化值.

> integer,dimension(2,2)::a = reshape((/1,2,3,4/),(/2,2/))

### 对数组使用内置函数

> allocated(array)   判断可分配数组的分配状态
>
> lbound(array)   给出数组array最小的索引指标
>
> ubound(array)  给出数组array最大的索引指标(integer,dimension(2,5)::a,lbound(a) –> 1,1 ;  ubound(a)–>(2,5)
>
> shape(array)	给出array的维数    shape(a)--->  2,5
>
> size(array)	给出array数组的大小   size(a)--->10
>
> all(array)	如果数组array中的所有元素都为真,逻辑函数返回.True.
>
> any(array)	如果数组array中的任意元素为真,逻辑函数返回.True.
>
> count(array)	统计array中真元素的个数
>
> dot_product(a1,a2)	计算两个大小相同的向量的点积
>
> matmul(a1,a2)	对两个一致的矩阵执行矩阵乘法
>
> maxloc(array,mask)	返回mask为真对应的array中的元素的最大值的位置,结果是带有一个元素的一维数组
>
> maxval(array,mask)	返回mask为真对应的array中元素的最大值(mask可以省略)
>
> **对于max*形式的函数,同样有min的函数来对应的求最小值的一些结果**
>
> sum(array,mask)	计算mask为真array中对应元素的和
>
> transpose(array)	返回array的转置矩阵

### 可分配数组

当一个问题所需要的数组大小不明确时,可以通过动态数组设置,在程序运行过程中,根据需求来动态的声明所需要的数组大小,其一可以避免内存的浪费,其二也可以避免由于数组大小补足而导致的程序错误.

> real,allocatable,dimension(:,:)::a	!声明一个大小为确定的二维数组
>
> allocate(a(10,10),stat=status)	!给数组一个确定的大小,status中保存执行该操作返回的一个结果,当其结果为0时,说明该操作成功的执行了,否则说明数组空间申请失败.
>
> **数组没有进行确定的大小分配时,是不可以使用它的**
>
> allocated(a)	!同样,在使用之前,可以通过来判断这个可变大小的数组,是否已经成功申请到的确定大小的内存空间,这样的判断可以避免程序中为申请空间而使用数组所产生一些错误.
>
> deallocate(a,stat=status)	!在使用完动态数组之后,要对申请的内存空间进行释放,避免浪费内存空间

> allocate(a1,source=a2,stat=status,errmsg=err_msg)	!这种形式的空间数组内存空间申请,表明a1和a2具有相同的大小,同样status保存申请成功与否的标志,err_msg中保存申请时遇到的错误信息,可以通过它来很好的判断程序的错误问题.

## 8.过程的附加特性

> real,save::sums	! save属性出现在类型声明语句中,在调用过程中任何一个用save属性定义的局部变量都会被保存,而不改变.
>
> program ex01
> implicit none
> integer i1,i2
> do i1 = 1,5
> 	call sub1()
> end do
> stop
> end
>
> subroutine sub1()
> implicit none
> integer,save::x1 = 0
> x1 = x1 + 1
> write(*,*)x1
> end subroutine sub1
>
> 在这个示例程序中,子过程sub1中x1使用的save属性,所以每一次循环调用sub1的时候,x1作为局部变量它的值是保存的,所以结果是  1,2,3,4,5

save属性,在通过module进行数据共享的时候,是一个必选的属性,这样就可以保证了,在一个地方对数据的修改,在另外的地方使用module中的共享数据时可以知道该数据的变化.

### 过程中的可分配数组

在一个子过程中,如果给一个可分配的数组使用save属性,那么在这个数组就只会在第一次使用的时候利用allocate分配数组大小,然后进行相关的计算,之后如果重复的使用该子过程,则不会重新进行内存的重新分配.相反的,如果没有这个save属性的利用,则在每次重复的利用该子过程的时候,程序都会重新进行动态内存的分配,然后进行计算,在计算结束之后,会将这个分配的空间释放掉.

## 9. 字符变量的特

字符变量之间是可以比较大小的，这种大小的比较是通过他们的ASCII码值进行的（通常都是选用ASCII码编码，但是如果用别的编码方式，那么比较的结果可能不相同）

> LLT 	字符串小于
>
> LLE	字符串小于等于
>
> LGT	字符串大于
>
> LGE	字符串大于等于
>
> 上面的这四个字符串函数，都是依照ASCII码值进行比较的，这样的程序可移植性强	

### 常用内置字符串函数

> achar(val)	给定一个val的整数值，返回它对应的ASCII字符
>
> char(val)	给定一个val的整数值，返回本机器使用编码库所对应的字符(这个结果依赖于机器所选编码)
>
> Iachar(string)  Ichar(string) 是上面连个函数的逆过程，由字符得到对应的整数值
>
> len(str)	返回字符串长度(它返回的是这个字符串在定义时候的长度)
>
> len_trim(str)	返回字符串的长度，不包括字符串尾部的空格(返回的是实际的非空格字符)
>
> trim(str)	删除掉字符串尾部的空格后，返回其结果  
>
> index(str1,str2,back)	用来查找str1中出现str2的索引位置

如果在一个子程序中要处理字符串，一般情况下字符串的长度是比较灵活的，所以可以声明未定长度的字符串

> character(len = *),dimension(n)::string1

通过这样的声明，就可以接受任意长度的字符串来进行处理，如果在这个子过程中需要使用字符串对应的长度，那么就可以通过len函数来获取在运行过程中传入的字符串的长度

> 可以这么做是因为，在声明的时候，既然是作为子过程的形式参数，那么它自然不过分配内存空间，当在调用的时候，有参数传递给子过程，我们知道Fortran中的参数传递是按地址进行的，所以就可以达到动态调整的目的

使用内存缓冲区，可以将字符形式的数据，变成数的形式

> character(len = 5)::input = '123.4'
>
> real::number
>
> read(input,*)number
>
> 这样就可以将一个字符类型的量成功的转换成数据形式

## 10.附加的内置数据类型

在通常的计算过程中，我们一般使用默认的数据类型(8字节)，这种默认的数据类型是存在一定的精度的，当我们向增加计算的精度时，可能默认的数据类型不能满足这种精度需求，所以我们就可以通过改变声明时候的数据类型来达到我们提高计算精度的目的。**我们可以使用模块的方式，来确定不同精度所使用的字节数，将它作为参数作为我们数据定义时候长度设置的变量，这样的话，如果我们想要提高我们的计算精度，只需要将这个模块中的参量进行修改，就可以达到目的**

> integer,parameter::single = 8
>
> integer,parameter::double = 16
>
> real(kind = single)::val1
>
> real(kind = double),dimension(20)::matrix

**kind函数可以用来判断一个量在计算机上存储时所占用的字节数**

### 选择与机器无关的数据精度方式

selected_real_kind(p = precision,r = range)，执行该函数的时候，返回适合或者超过指定取值范围和精度的实数类型的最小类别的类别号。函数中的参数precision是所需精度的十进制数，range是以10的幂形式所表示的所需要的指数范围。通俗点讲，我们以科学计数法的角度来表示数据，precision就是**有效数字**，range是**幂指数的大小**

**可以同过给定一个计算过程中，你认为最大或者最小的数的科学计数法的形式，那么该函数就可以返回一个正整数，告诉你该在该计算机上需要用多少字节来才能表示这个数**

selected_int_kind(r = range)，该函数是对整数型变量操作，range的含义和上面一致

selected_char_kind(name)	该函数作用于字符串上，给它输入一个字符串，判断需要使用何种类型参数来合理的存储这个字符串形式

通过这个函数确定后的数据类型的精度，使得程序具有很大的可移植性。

> range(x)	该函数会判断x的数据类型是单精度还是双精度，然后返回在这个精度下，数据以科学计数法表示时，它幂指数的最大值是多少
>
> precision(x)	该函数返回数据x的有效数字
>
> dble(x)	无论x是何种精度的数据，都将其进行类型变换，使其成为双精度的数据类型

> 复数是由两部分组成：实部和虚部，所以如果声明一个复数是8个字节的，那么其实它的由16个字节来存储的，分别要保存它的实部和虚部
>
> cabs(val)	返回一个复数量val的辐值，也就是复数的模长
>
> cmplx(a,b,kind)	通过a和b来构建一个需要的复数变量
>
> int(val)	将复数量val的实部转换成整数

## 11.过程和模块的高级特性

### 递归过程

在Fortran中，通常一个函数是不可以调用自己本身的，然而一些实际的问题如果用递归调用的方法则会更加简单，比如计算阶乘。如果要让一个函数可以调用自己本身，那么需要在定义这个函数的时候，明确这是一个递归函数，那么Fortran编译器则会以直接或者间接的方式来调用函数自身。

> recursive subroutine fact(n,result)
>
> implicit none
>
> integer,intent(in)::n
>
> integer,intent(out)::result
>
> integer temp
>
> if(n>1)then
>
> ​	call fact(n-1,temp)
>
> ​	result = n*temp
>
> else
>
> ​	result  = 1
>
> end if
>
> end subroutine fact

如上所示定义了一个递归的子过程，从而实现了求阶乘。在Fortran中，我们知道，子过程是不会由返回值的，如果向通过一个子过程得到计算的结果，都是通过子过程的形参列表来进行的。然而函数不一样，函数是可以有返回值的，所以上述求阶乘同样可以通过定义递归的函数来实现，只不过在实现递归函数的时候有一些问题。

如果函数要调用自己，当设置返回值的时候，在这个递归的过程中，函数名无法避免的要出现在等号(=)的左侧，这样的话会导致函数名，既当作函数名，又当做变量，会产生混乱，所以，我们可以通过向某个特定的形参指定为返回值，这样就避免了函数名作用冲突的矛盾。

> recursive function fact(n) result(answer)
>
> implicit none
>
> integer,intent(in)::n
>
> integer::answer
>
> if(n>1)then
>
> ​	answer = n*fact(n-1)
>
> else
>
> ​	answer = 1
>
> end if
>
> end function fact

其实本质上，使用递归函数的方法只不过是借用了子过程的方式，通过一个形式参数来当做函数的返回值。

## 12.过程接口和接口块

通常情况下，给一个子过程创建一个接口的方式，是将子过程写在一个模块之中，但是这样的方法并不是很方便。当不可能将子过程放在模块中的时候，可以在调用这个子过程的单元中定义接口快，在定义的这个接口块之中，它包含了所要调用的子过程的所有接口特征，Fortran编译器则会根据所提供的接口信息执行一致性检查

> interface
>
> ​	interface_boday1
>
> ​	interface_boday2
>
> end interface

接口块的一般结构如上所示，其中每个接口块中都是由外部过程的初始化subroutine和function语句，过程参数相关的特定类型声明语句。

> program interface_example
>
> implicit none
>
> interface 
>
> ​	subroutine sort (a,n)
>
> ​	implicit none
>
> ​	real,dimension(:),intent(inout)::a
>
> ​	integer,intent(in)::n
>
> ​	end subroutine sort
>
> end interface

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