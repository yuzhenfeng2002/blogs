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

## 模型设置

- 状态 $\omega \in \Omega$
- 信息接收者
  - 行动 $a \in A$
  - 效用 $u(a, \omega)$
  - 对状态的先验分布 $\mu_0$
  - 对状态的后验分布 $\mu_s$
  - 后验分布的分布 $\tau$
  - 最优行动 $a(\mu) = \arg\max \mathbb{E}_\mu [u(a, \omega)]$
    - 如果期望相等，则会选择使信息发送者效用更高的行动
- 信息发送者
  - 效用 $v(a, \omega)$
  - 信号
    - 策略 $\pi(\cdot|\omega)$
    - 实现 $s \in S$
    - （承诺假设）信息发送者的信号被真实地传递给接收者；发送者可以隐藏部分或全部信息，只要他发送的信号是可验证的。

根据贝叶斯公式，后验分布
$$
\mu_s(\omega) = \frac{\pi(s|\omega) \mu_0(\omega)}{\sum_{\omega'\in \Omega} \pi(s|\omega') \mu_0(\omega')}
$$
后验分布的分布
$$
\tau(\mu) = \sum_{s\in S} \mathbb{I}(\mu_s = \mu) \sum_{\omega'\in \Omega} \pi(s|\omega') \mu_0(\omega')
$$
后验分布的分布的均值等于先验分布
$$
\begin{aligned}
\sum_{\mu}\mu(\omega) \tau(\mu) = \sum_{\mu} \sum_{s\in S} \mathbb{I}(\mu_s = \mu) \mu(\omega) \sum_{\omega'\in \Omega} \pi(s|\omega') \mu_0(\omega')
\\ = \sum_{\mu} \sum_{s\in S} \mathbb{I}(\mu_s = \mu) \pi(s|\omega) \mu_0(\omega)
\\ = \mu_0(\omega) \sum_{\mu} \mathbb{E}[\mathbb{I}(\mu_{\tilde{s}} = \mu)|\omega]
\\ = \mu_0(\omega) \sum_{\mu} \operatorname{Pr}\{\mu_{\tilde{s}} = \mu|\omega\} = \mu_0(\omega)
\end{aligned}
$$

### 简化问题（专注于后验分布的分布）

后验分布的分布是由$\pi$诱导产生的，如果不同的$\pi$诱导产生了同样的$\tau$，那么这不同的$\pi$在给定状态的条件下对信息发送者的效用是相同的。我们还需要知道，什么样的$\tau$可以被诱导，以及$\tau$对应的效用是多少。

我们先考察$\tau$对应的效用。给定后验分布为$\mu$的条件下，发送者的效用为
$$
v(\mu) = \mathbb{E}_{\mu_0}[v(a(\mu), \omega)|\mu] = \frac{\sum_{w\in \Omega} \sum_{s\in S} \mathbb{I}(\mu_s=\mu) v(a(\mu),\omega) \pi(s|\omega)\mu_0(\omega)}{\sum_{w\in \Omega} \sum_{s\in S} \mathbb{I}(\mu_s=\mu) \pi(s|\omega)\mu_0(\omega)}
$$
其中，分母即为$\tau(\mu)$，分子为
$$
\sum_{w\in \Omega} v(a(\mu),\omega) \mu(\omega) \sum_{s\in S} \mathbb{I}(\mu_s=\mu) \sum_{\omega'\in \Omega} \pi(s|\omega') \mu_0(\omega')
= \sum_{w\in \Omega} v(a(\mu),\omega) \mu(\omega) \tau(\mu)
= \tau(\mu) \mathbb{E}_{\mu}[v(a(\mu), \omega)]
$$
因此
$$
v(\mu) = \mathbb{E}_{\mu_0}[v(a(\mu), \omega)|\mu] = \mathbb{E}_{\mu}[v(a(\mu), \omega)]
$$
进而，后验分布的分布$\tau$产生的平均效用为
$$
\mathbb{E}_{\tau}[v(\mu)] = \sum_{\mu} v(\mu) \tau(\mu)
$$

接着，为了引出定理 1，我们定义
> 信号$\pi$的值
> $$
  \sum_{w\in \Omega} \sum_{s\in S} \pi(s|\omega) v(a(\mu_s), \omega) \mu_0(\omega)
  $$

> 贝叶斯可说服的后验分布的分布 $\tau$
> $$
  \sum_{\mu} \mu \tau(\mu) = \mu_0
  $$

> 直接的信号
> $$
  S\subseteq A,~ a(\mu_s) = s
  $$

> 定理 1
> 
> 以下三者等价
> 1. 存在信号策略 $\pi$，使得$\pi$的值为 $v^*$；
> 2. 存在**直接的**信号策略 $\pi$，使得$\pi$的值为 $v^*$；
> 3. 存在**贝叶斯可说服的**后验分布的分布 $\tau$，使得 $\mathbb{E}_{\tau}[v(\mu)] = v^*$。

有了定理 1，我们一方面可以只关注直接的信号，另一方面，研究最优信号就只要解决
$$
\max \quad \mathbb{E}_{\tau}[v(\mu)] \text{ s.t. } \sum_{\mu} \mu \tau(\mu) = \mu_0
$$

> 定理 1 证明
> 
> 根据定义有 (2)->(1)，根据上面的结论“后验分布的分布的均值等于先验分布”有 (1)->(3)。
>
> 我们先证明 (1)->(2)，这是通过根据原信号构造了新的“直接的”信号。如果接收者基于原信号的行动是$a$，那么新的信号就直接“建议”$a$。[数学上如何说明？](#)
>
> 我们再证明 (3)->(1)，这也是通过构造证明的。令$\pi(s|\omega) = \mu_s(\omega) \tau(\mu_s) / \mu_0(\omega)$，则
> $$
  \begin{aligned}
    \sum_{w\in \Omega} \sum_{s\in S} \pi(s|\omega) v(a(\mu_s), \omega) \mu_0(\omega)\\ 
    = \sum_{w\in \Omega}\sum_{s\in S} v(a(\mu_s), \omega) \mu_s(\omega) \tau(\mu_s) \\
    = \sum_{s\in S} \tau(\mu_s) v(\mu_s) = \sum_{\mu} \tau(\mu) v(\mu) = v^*
  \end{aligned}
  $$
> 贝叶斯可说服保证了
> $$
  \sum_{s\in S} \pi(s|\omega) = \sum_{s\in S} \mu_s(\omega) \tau(\mu_s) / \mu_0(\omega) = 1
  $$

下面，我们只需要在$\mu$的空间里研究$v(\mu)$。对于任意$(\mu^*, v(\mu^*))$，只要它在$v(\mu)$下境图的凸包里，就能通过一些$(\mu, v(\mu))$凸组合构造，因此我们也只需要研究这个凸包。定义
> $$
  V(\mu) = \sup \{z|(\mu, -z)\in \operatorname{conv}(\operatorname{epi}(-v(\mu)))\}
  $$
根据上述分析，我们有推论
> 最优信号的值为
>
> $$
  V(\mu_0) = \max \quad \mathbb{E}_{\tau}[v(\mu)] \text{ s.t. } \sum_{\mu} \mu \tau(\mu) = \mu_0
  $$

## 什么时候能从说服中获利？