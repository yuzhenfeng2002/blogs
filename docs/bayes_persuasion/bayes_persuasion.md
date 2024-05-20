# Bayesian Persuasion

Kamenica, E., & Gentzkow, M. (2011). Bayesian persuasion. American Economic Review, 101(6), 2590-2615.

## 引例

有公诉人与法官。公诉人想要说服法官判被告有罪，法官想要进行正确的判决。公诉人和法官有先验知识被告有罪的概率为$p<\frac{1}{2}$；公诉人进行了调查，会依概率$s, t$给法官如下信息$I$：

|        |      | 被告 |有罪(g)|无罪(i)|
|:------:|:----:|:----:|:----:|:----:|
|公诉人(I)| 有罪 |      | $s$ | $1-t$ |
|        | 无罪 |      | $1-s$ | $t$ |

法官根据信息$I$，有了后验概率
$$
P_{g|g} = P\{g|I=g\} = \frac{P\{I=g|g\}P\{g\}}{P\{I=g\}} = \frac{ps}{ps+(1-p)(1-t)}
$$
$$
P_{i|i} = P\{i|I=i\} = \frac{P\{I=i|i\}P\{i\}}{P\{I=i\}} = \frac{(1-p)t}{p(1-s)+(1-p)t}$$

由于法官想要进行正确的判决，因此，当$I=g$时，他会选择$P_{g|g}$和$P_{i|g} = 1 - P_{g|g}$中的较大者，$I=i$时同理。他的平均效用是
$$
p(s \mathbb{I}(P_{g|g} \geq P_{i|g}) + (1-s)\mathbb{I}(P_{g|i} \geq P_{i|i})) + (1-p)((1-t) \mathbb{I}(P_{g|g} < P_{i|g}) + t\mathbb{I}(P_{g|i} < P_{i|i}))
$$

公诉人想要说服法官判被告有罪，他的平均效用是
$$
\begin{aligned}
p(s \mathbb{I}(P_{g|g} \geq P_{i|g}) + (1-s)\mathbb{I}(P_{g|i} \geq P_{i|i})) + (1-p)((1-t) \mathbb{I}(P_{g|g} \geq P_{i|g}) + t\mathbb{I}(P_{g|i} \geq P_{i|i}))\\
= (ps+(1-p)(1-t)) \mathbb{I}(ps\geq(1-p)(1-t)) + (p(1-s)+(1-p)t) \mathbb{I}(p(1-s)\geq(1-p)t)\\
= 2ps \mathbb{I}(ps\geq(1-p)(1-t)) + 2p(1-s) \mathbb{I}(p(1-s)\geq(1-p)t)
\end{aligned}
$$
以上两个示性函数里的区域是不相交的，最优值分别可以取在
$$
(s,t)=(1, 1-\frac{p}{1-p}), \quad (s,t)=(0, \frac{p}{1-p})
$$
由于公诉人不应该隐瞒犯罪事实，所以我们取前者，公诉人的最优效用为$2p$。

可以看到，如果公诉人完全不提供证据，法官会判所有人无罪，由于$p<\frac{1}{2}$；如果公诉人提供完全的信息，法官会判$p$比例的人有罪；如果公诉人进行如上的披露，法官会判$2p$比例的人有罪，公诉人的效用得到了提升。