# 用经济学解释内卷（二）：你卷赢了吗？

书接上回：
> 假设你面前有两条路，一条是努力，一条是躺平。由于资源的限制，努力的你不一定会收获好的结果，你会选择什么？

与上回不同，这次的你可以选择努力的程度，和对手竞争有限的资源。你的努力，会为你带来什么？

> **请填写关联词**
> 
> __ __ 累死自己，__ __ 卷死别人。

![](network.png)

我们继续将这个问题建模为一个博弈问题，如下图所示。假设有 $N$ 个人，每个人可以选择努力的程度 $e$。躺平时，$e = 0$，收获 $u_0$；努力时，$e > 0$，但只有努力程度排名前 $n$ 的人才能收获 $u_1$；剩余的人只能和躺平的人一样，收获 $u_0$。记排名第 $n$ 的人的努力程度为 $e_{(n)}$，则一个人的效用可以计算为：
$$
  u(e) = \begin{cases}
    u_0 - e,  & 0 \leq e < e_{(n)}\\
    u_1 - e,  & e \geq e_{(n)}\\
  \end{cases}
$$

资源充足时，所有人都会小小努力一点点，然后收获 $u_1$。因此，我们主要讨论资源有限时的情况，这也是如今的现实。

在均衡时，如果某人付出的努力 $e > e_{(n)}$，则他可以减少努力获得更多的效用；如果某人付出的努力 $e < e_{(n)}$，则他可以更加努力或直接躺平获得更多的效用。因此，在均衡时，一个人付出的努力要么是 $0$，要么是 $e_{(n)}$。如果 $e_{(n)} < u_1 - u_0$，就会有人比第 $n$ 个人多努力那么一点点，从而获得更多的效用，这时 $e_{(n)}$ 增大；如果 $e_{(n)} > u_1 - u_0$，那么就有人觉得得不偿失，放弃努力，这时 $e_{(n)}$ 减小。

因此，均衡时，平均会有 $n$ 个人选择付出 $u_1 - u_0$ 的努力，剩余的人会选择躺平。但，我们会发现，**努力后的效用和直接躺平是一样的**。这也是如今大家的现实困境：努力了，但似乎收益甚微。

那么，所有所谓的上岸，所谓高考、考研、考公，当你到达渴望的彼岸，你真的卷赢了吗？或许我们应该乐观一点，既然无论做什么选择都很难从这个世界套利，那么不如 Just do it，Just follow my heart。