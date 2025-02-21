---
title: 科研日志模板(Latex)
tags: Topology Latex
layout: article
license: true
toc: true
key: a20220909
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
平时在阅读文献或者组会讨论的时候总会有一些需要记录的内容，这里就修改了别人的一个模板，整理了一个适合自己的记录科研log的[Latex模板](https://github.com/yxli8023/Research-Log)。
{:.info}
<!--more-->

# Research-Log
**这是用来整理平时科研过程中自己的想法，和别人合作时各方的讨论结果以及问题整理的一个模版。平时在开组会活着与同学讨论过程中总会有一些想法，顺便可以整理记录，时间久了应该也能从中找出一些不错的点子。**

我这里的模版是基于[project-logbook](https://github.com/apalha/project-logbook)进行了一些改造。它是一个英文版本，就是用来写科研日志的，我在这里进行了中文化，并对其中的一些环境进行了修改，对高亮色块的颜色配置进行了修改。

# Meeting
首先是**Meeting**环境
```latex
\begin{Meeting}{会议日期}{main}%最后一个参数是作者id，用来标记这次会议是谁的想法
这里是一次会议讨论内容的记录，主导者第一个作者，用mian来作为他的id。在这次讨论中有一个重要任务先要完成
\hightodo{日期}{main}{这里是目前最优先要做的事情，这些事情都会被加入到最后的todo list 中 }
\end{Meeting}
```
这个环境块主要是用来整理一下会议或者讨论的记录，其中可以标记出会议的主导者。同时也能利用
```latex
\hightodo{日期}{main}{这里是目前最优先要做的事情，这些事情都会被加入到最后的todo list 中 }
```
高亮出讨论中最重要，优先级最高的需要解决的问题。

# Note(Think)
这个环境块和**Meeting**没有太大本质上的区别，就是进行了颜色上的区分，可以给这个**Note**环境一个标题，用来提醒自己这个Note要说明什么问题。同样还有一个**Think**环境块，也只是修改了一下颜色而已。

这里还有一个**lotodo**，用来区别前面的**hightodo**，二者的颜色是不同的，代表了需要完成的事情的优先级是不同的，相当于是要高亮一下提醒自己需要做的事情，
```latex
\lowtodo{日期}{abc}{这里是目前次优先要做的事情 ，这些事情都会被加入到最后的todo list 中}
```
![png](/assets/images/logo/demo-1.png)

最后在文档最后面会生成一个**Todo List**，这里面就包含了前面使用**hightodo**和**lowtodo**来提醒自己需要做的事情，方便自己在整理科研日志的时候来提醒自己，同时也可以将一些想法高亮出来，说不定什么时候就可以有一些突破性的进展。
![png](/assets/images/logo/demo-2.png)

# 未待续完
暂时这个模版满足了我的基本需求，还没有其他的想法可以加入到这个模版中，后面如果有新东西，我会继续更新。如果你在科研进程中有什么好的整理日志的方式，可以联系我，如果在我能力范围内，我会努力将这个功能加入到当前的模版中。




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