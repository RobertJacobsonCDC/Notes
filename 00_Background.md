# 0.0: Basic Definitions from College Probability and Statistics 101


- **Probability Density Function (PDF)**: A function describing the relative likelihood of a continuous random variable taking a specific value.

For continuous distributions, this definition isn't good enough, because $`P(X=c)=0`$ for all numbers $`c`$! Instead the PDF describes the probability that an event will occur in some arbitrary interval. If $`f_X(x)`$ is the PDF for the distribution of a random variable $`X`$, then

```math
P(a \lt X\lt b)=\int_a^bf_X(x) dx
```

**Observation:** *If you know everything about a probability distribution on intervals, then you know the distribution itself.*

**Observation:** For continuous probability distributions, $`P(a\leq X\lt b)=P(a\lt X\leq b)=P(a\leq X \leq b)=P(a\lt X\lt b)`$, because the probability of $`X=c`$ is exactly zero; the probability is spread out over an infinite number of values $`c`$.



- **Cumulative Distribution Function (CDF)**: A function giving the probability that a random variable is less than or equal to a specific value. For a continuous random variable with PDF $`f_X`$, $`F_X(x) := P(X \leq x) = \int_{-\infty}^x f_X(w) dw`$.

The CDF tells you probabilities on "half intervals" $`(-\infty, x]`$. But it turns out that if you know probabilities on "half intervals," you know everything about the distribution itself. First, notice we automatically have probabilities on any interval $`(a, b)`$:

```math
P(a\lt X\lt b)=\int_a^bf_X(x) dx=\int_{-\infty}^b f_X(x) dx-\int_{-\infty}^a f_X(x) dx=F_X(b)-F_X(a)
```

**Observation:** We have $`F_X(x)- F_X(c)=\int_c^xf_X(t) dt`$. The Fundamental Theorem of Calculus says that $`\frac{d}{dx}\left[\int_c^xf_X(t) dt\right]=f_X(x)`$. So,

```math
\frac{d}{dx}\left[F_X(x)- F_X(c)\right]=\frac{d}{dx}\left[\int_c^xf_X(t) dt\right]=f_X(x).
```

The RHS is independent of $`c`$, so the LHS is, too. (Let $`c\to -\infty`$.) This shows:

**Theorem:** The PDF $`f_X(t)`$ and CDF $`F_X(t)`$ are related by $`\frac{d}{dt}\left[F_X(t)\right]=f_X(t)`$.



- **Survival Function**: A function giving the probability that a random variable *exceeds* a specific value, $`S_X(x) = P(X > x) = 1 - F_X(x)`$.
- **Hazard Function**: A function describing the instantaneous rate at which an event occurs at time $`t`$, given that the event has not occurred up to time $`t`$:

```math
h_X(t):=\lim_{\Delta t\to0}\frac{P(t\leq X \leq t+\Delta t \mid X\geq t)}{\Delta t}.
```

**Theorem:** $`h_X(t)=\frac{f_X(t)}{S_X(t)}`$.

**Proof:**

$`h_X(t):=\lim_{\Delta t\to0}\frac{P(t\leq X \leq t+\Delta t \mid X\geq t)}{\Delta t}`$

$`\hspace{2cm}=\lim_{\Delta t\to0}\frac{P(t\leq X \leq t+\Delta t)}{P(X\geq t) \Delta t}`$

$`\hspace{2cm}=\frac{\lim_{\Delta t\to0}\frac{P(t\leq X \leq t+\Delta t)}{\Delta t}}{P(X\geq t)}`$

$`\hspace{2cm}=\frac{\frac{dF_X}{dt}(t)}{S_X(t)}`$

$`\hspace{2cm}=\frac{f_X(t)}{S_X(t)}\hspace{5cm}\Box`$


With $`S_X(x) = P(X > x) = 1 -  P(X \leq x) = 1 - F_X(x)`$, then we also have 
```math
P(a\lt X\lt b)= F_X(b) - F_X(a) = S_X(a) - S_X(b)
```


# 0.1: Background: Poisson Processes and the Exponential Distribution

We will work up to the proof of the following theorem.

> **Theorem:** The only continuous memoryless probability distribution is the exponential distribution.

## Poisson and Exponential Distributions

First recall some definitions.


- **Homogeneous Poisson Process**: A stochastic counting process that models events occurring
independently at a constant average rate $`\lambda`$.

- **Poisson Distribution**: A *discrete* probability distribution giving the probability of $`k`$ events occurring by time $`t`$ , with probability mass function (PMF; for discrete distributions, analogous to the PDF for continuous distributions) of

```math
p_{N(t)}(k) = P(N(t) = k) = \frac{(\lambda t)^k e^{-\lambda t}}{k!}, \qquad k=0,1,2,\ldots
```

- **Exponential Distribution**: A continuous probability distribution describing the time between events in a Poisson process, with PDF $`f(x) = \lambda e^{-\lambda x}`$ for $`x \geq 0`$.  If a random variable $`X`$ has exponential distribution with rate $`\lambda`$, we write $`X \sim \text{Exp}(\lambda)`$.


## Memorylessness

 - The distribution of a positive continuously distributed random variable $`X`$ is **memoryless** if $`P(X>t) = P(X > s+t \mid X > s)`$ for all $`s,  t\geq 0`$.

In words, the probability that an event occurs after a time $`t`$ (from the start) is the same as the probability that it occurs after an additional time $`t`$ from a time $`s`$ given that it has not occurred by time $`s`$. The time we have to wait for an event is independent of the time we have already waited.

**Lemma:** The distribution of a positive continuously distributed random variable $`X`$ is memoryless iff $`P(0\lt  X \lt  t) = P(s\lt  X \lt  t+s \mid X > s)`$.

**Proof:**

$`\hspace{0.5cm}P(X>t) = P(X > s+t \mid X > s)`$

$`\to 1-P(X>t) = 1-P(X > s+t \mid X > s)`$


$`\to P(0\lt  X \lt  t) = P(s\lt  X \lt  t+s \mid X > s)`$

Conversely, suppose that $`P(0\lt X\lt t)=P(s\lt X\lt s+t\mid X>s)`$. Since $`X`$ is positive and continuously distributed, endpoint probabilities are zero. Thus,

```math
\begin{aligned}
P(X>t)
&=1-P(0\lt X\lt t)\\
&=1-P(s\lt X\lt s+t\mid X>s)\\
&=P(X>s+t\mid X>s).
\end{aligned}
```

Therefore $`X`$ is memoryless. $`\hspace{5cm}\Box`$



## Main Theorem

We are ready to prove the main theorem of this section.

**Theorem:** The only continuous memoryless probability distribution is the exponential distribution.

**Proof:** For every $`\Delta t>0`$, memorylessness gives

```math
P(t<X\leq t+\Delta t\mid X>t)=P(0<X\leq\Delta t).
```

Therefore,

```math
h(t)=\lim_{\Delta t\to0}\frac{P(t<X\leq t+\Delta t\mid X>t)}{\Delta t}
=\lim_{\Delta t\to0}\frac{P(0<X\leq\Delta t)}{\Delta t}=\lambda,
```

which is independent of $t$. Thus the hazard is a constant $`\lambda>0`$. Also, since $`S(t)=1-F(t)`$, we have $`\frac{dS}{dt} = -\frac{dF}{dt}`$. It follows that

```math
\lambda = h(t) = \frac{\frac{dF}{dt}}{S(t)} = \frac{-\frac{dS}{dt}}{S(t)} \longleftrightarrow  -\lambda S(t) = \frac{dS}{dt} \longleftrightarrow S(t)=e^{-\lambda t}+ c.
```

The initial condition $`S(0)=P(X\geq 0) = 1`$ gives $`c=0`$. Since $`\lambda =h(t) = \frac{f(t)}{S(t)}`$, it follows that the PDF function $`f(t) = \lambda e^{-\lambda t} \quad (t\geq 0). \hspace{5cm}\Box`$

**Alternative Proof:** The memoryless condition says

$`P(X>t) = P(X > s+t \mid X> s) = \frac{P(X > s+t)}{P(X > s)}`$

$`\to P(X > s+t) = P(X>t)\cdot P(X>s)`$

To apply this repeatedly, let $`s=(a-1)t`$. Then

```math
P(X>at)=P(X>(a-1)t+t)=P(X>(a-1)t)P(X>t).
```

For example,

```math
P(X>3t)=P(X>2t)P(X>t)=P(X>t)^2P(X>t)=P(X>t)^3.
```

By induction, $`P(X>at)=P(X>t)^a`$ for every $`a\in\mathbb{N}`$.

For $`b\in\mathbb{N}`$, apply the integer result with $`t/b`$ in place of $`t`$:

```math
P(X>t)=P\left(X>b\frac{t}{b}\right)=P\left(X>\frac{t}{b}\right)^b.
```

Thus $`P(X>t/b)=P(X>t)^{1/b}`$. Combining this with the integer result, if $`q=a/b\in\mathbb{Q}_{\geq0}`$, then

```math
P(X>qt)=P\left(X>a\frac{t}{b}\right)=P\left(X>\frac{t}{b}\right)^a=P(X>t)^q.
```

Finally, the survival function $`S(t)=P(X>t)`$ is continuous because $`X`$ has a continuous distribution. Since rational numbers can approximate every real number, this equality extends to $`P(X>rt)=P(X>t)^r`$ for all $`r\geq0`$.

Taking $`t=1`$, we have $`P(X>r)=P(X>1)^r`$. Therefore $`P(X>t)=e^{-\lambda t}`$, where $`\lambda=-\ln(P(X>1))=-\ln(S(1))`$. Thus, we have shown that $`S(t)=e^{-\lambda t}`$ for some $`\lambda>0`$. The CDF is thus $`F(t)=1-S(t)=1-e^{-\lambda t}`$. The PDF is $`f(t)=\frac{dF}{dt}(t)=\lambda e^{-\lambda t}.\hspace{5cm}\Box`$


# 0.2: Time Scaling

 #### A Poisson process is **inhomogeneous** if the rate $`\lambda`$ is not constant over time.

Let $`\lambda=\lambda(t)`$ be a non-negative integrable rate function, and let $`T`$ be a waiting-time random variable with hazard $`\lambda(t)`$. For finite $`t`$, its density is $`f_T(t)=\lambda(t)e^{-\int_0^t\lambda(s)ds}`$, and $`F_T(t)=P(T\leq t)=1-e^{-\int_0^t\lambda(s)ds}`$. This distribution is often called the waiting-time or first-arrival-time distribution.  Let $`\Lambda(t)`$ be the cumulative rate function defined by

$`\displaystyle\hspace{2cm}  \Lambda(t) := \int_0^t \lambda(s)  ds`$

Then we can write

$`\displaystyle\hspace{2cm} f_T(t) = \lambda(t) e^{-\Lambda(t)}`$, and

$`\displaystyle\hspace{2cm} F_T(t) = 1 - e^{-\Lambda(t)}`$.


If $`\lambda(t)>0`$, then $`\Lambda`$ is strictly increasing and has an ordinary inverse on its range. Note, that when $`\lambda`$ is constant and positive, $`\Lambda(t) = \lambda t`$ and the waiting-time distribution simplifies to the usual exponential distribution.

If $`\lambda`$ is zero on an interval, $`\Lambda`$ is flat there; in that case use its generalized inverse $`\Lambda^{-1}(y):=\inf\{t\geq0:\Lambda(t)\geq y\}`$.

If $`\Lambda(\infty)<\infty`$, then there is a nonzero probability $`e^{-\Lambda(\infty)}`$ that no event ever occurs; equivalently, take $`T=\infty`$ on that event. This is the case for a finite infectious period.



**Theorem:** Let $`\lambda=\lambda(t)>0`$ be a rate function with $`\Lambda(\infty)=\infty`$, and let $`d=\Lambda^{-1}`$. If $`X\sim\operatorname{Exp}(1)`$, then $`d(X)`$ has hazard function $`\lambda(t)`$.

Equivalently, if $`T`$ has hazard $`\lambda(t)`$, then $`\Lambda(T)\sim\operatorname{Exp}(1)`$. The map $`\Lambda`$ sends a time to cumulative-rate space, and its inverse $`d`$ sends a cumulative-rate value back to its corresponding time. Thus, to sample $`T`$, sample $`X\sim\operatorname{Exp}(1)`$ in cumulative-rate space and set $`T=d(X)`$.

**Example:** Suppose $`\lambda(t)=0.7`$. Then $`\Lambda(t)=\int_0^t0.7\,ds=0.7t`$, so its inverse is $`d(y)=\frac{y}{0.7}`$. The theorem says to start with $`X\sim\operatorname{Exp}(1)`$ in cumulative-rate space and map it back to time by applying $`d`$. Thus

```math
T=d(X)=\frac{X}{0.7}.
```

The resulting waiting time has constant hazard $`0.7`$, so $`T\sim\operatorname{Exp}(0.7)`$.

This says that for a Poisson process, if the rate ($`\lambda=0.7`$) is smaller (than $`1`$ in this case), the wait time until the next event ($`\frac{1}{0.7}X`$) is larger.

**Proof:** We perform a change of variables $`Y=d(X)`$. By definition,

$`F_Y(y) := P(Y \leq y) = P(d(X) \leq y)  \overset{\star}{=} P(X \leq \Lambda(y)) =: F_X(\Lambda(y))`$.

Now we use the Chain Rule to differentiate:

$`f_Y(y) := \frac{d}{dy}\left[ F_Y(y) \right] = \frac{d}{dy}\left[ F_X(\Lambda(y)) \right] = f_X(\Lambda(y)) \cdot \frac{d}{dy}\left[ \Lambda(y) \right]`$.

By the Fundamental Theorem of Calculus, $`\frac{d}{dy}\left[ \Lambda(y) \right] = \lambda(y)`$. Thus, $`f_Y(y) = \lambda(y)e^{-\Lambda(y)}`$, the PDF of a waiting time with hazard $`\lambda(y)`$. This shows that $`Y=d(X)`$ has hazard $`\lambda`$. $`\hspace{3cm}\Box`$

The equality marked with a '$`\star`$' needs justification, because it's not true in general. If $`a \lt  b`$, what conditions on a function $`f`$ guarantee that $`f(a) \leq f(b)`$? If ($`\star`$) is to hold, $`\Lambda(t)`$ must satisfy those conditions.

**Answer:** $`f(a) \leq f(b)`$ clearly holds if $`f`$ is monotonically increasing (by definition). Since $`\lambda(t)>0`$, $`\Lambda(t) := \int_0^t \lambda(s)  ds`$ is in fact strictly increasing.


**Another proof:**

Let $`\lambda = \lambda(t)`$ be some non-negative real-valued rate function of a Poisson point process, i.e. the hazard function. When $`\lambda(t)`$ is not constant, the Poisson point process is inhomogeneous.

Let $c(t)$ be the cumulative rate function, s.t.

```math
c(t) = \int_0^t \lambda(s)ds
```
With $`\lambda(t)`$ a non-negative real-valued function, the cumulative rate function is monotonically non-decreasing. It has an ordinary inverse only when it is strictly increasing; otherwise, use the generalized inverse defined above.

$`\lambda(t)`$ as the hazard function also satisfies the relationship,

```math
\lambda(t) = \frac{f(t)}{S(t)} = \frac{1}{S(t)} \frac{dF}{dt} = \frac{-1}{S(t)}\frac{dS}{dt} = -\frac{d}{dt}\Big(\ln S(t)\Big)
```

Then

```math
ln(S(t)) = -\int_0^t \lambda(s)ds = -c(t)
```

Therefore, we can write the survival function $`S(t) = P(T > t)`$ as

```math
S(t) = e^{-c(t)}
```

For the following ordinary-inverse argument, assume $`c`$ is strictly increasing. Then it has an inverse $`d`$ such that $`d(c(t))=t`$ and $`c(d(t))=t`$. If $`c`$ has flat portions, the generalized inverse gives the corresponding sampling construction.

Therefore, $`c(T) > y \Leftrightarrow T > d(y)`$ since $c(t)$ and $d(t)$ uniquely map to each other.

Then

```math
P(c(T) > y) = P(T > d(y)) = S(d(y)) = e^{-c(d(y))} = e^{-y} = S^{*}(y)
```

where $`S^{\text{*}}(y)`$ is the survival function of events in the cumulative space.

This means that $`c(T) = y \sim \operatorname{Exp}(1)`$, i.e. inter event distances in the cumulative space can be sampled with an exponential distribution of rate equal to $1$.

**And another proof:**

A waiting time $`T`$ with hazard $`\lambda(t)`$ has survival function

```math
P(T > t) = S_T(t) = e^{-\Lambda(t)}
```

Then, for $y \ge 0$,

```math
P(\Lambda(T) > y) = P(T > \Lambda^{-1}(y)) = S_T(\Lambda^{-1}(y)) = e^{-\Lambda(\Lambda^{-1}(y))} = e^{-y}.
```

But that is exactly the survival function of $`\text{Exp}(1)`$, so $`\Lambda(T) \sim \text{Exp}(1)`$.

Applying $`\Lambda^{-1}`$ to both sides means $`T = \Lambda^{-1}(X), X \sim \text{Exp}(1)`$.

---



Now, what does this mean for us in terms of being able to sample the time to the next infection when the rate process is inhomogeneous?

Since we can sample the interevent distances in the cumulative space with $\sim \text{Exp}(1)$, then we can get the cumulative values $y_i$ by summing the sampled inter event values $\Delta y_i$, and invert each $y_i$ with $d(y)$ to get the time for each event $t_i$.

Algorithmically, this means:

1. Sample $\Delta y_i$ from $\operatorname{Exp}(1)$
2. Calculate $y_i$ as $y_i = \Sigma_{i' = 1}^i \Delta y_{i'} $
3. Invert $y_i$ to get $t_i$ through $t_i = d(y_i)$
4. Calculate the inter event times $\Delta t_i = t_i - t_{i-1}$

where $\Delta t_1 = t_1 - t_0 = t_1 - 0 = t_1$ since the time of the zeroeth event would be 0.
