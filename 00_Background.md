# 0.0: Basic Definitions from College Probability and Statistics 101


- **Probability Density Function (PDF)**: A function describing the relative likelihood of a continuous random variable taking a specific value.

For continuous distributions, this definition isn't good enough, because $P(X=c)=0$ for all numbers $c$! Instead the PDF describes the probability that an event will occur in some arbitrary interval. If $f_X(x)$ is the PDF for the distribution of an r.v. $X$, then

```math
P(a \lt X\lt b)=\int_a^bf_X(x) dx
```

**Observation:** *If you know everything about a probability distribution on intervals, then you know the distribution itself.*

**Observation:** For continuous probability distributions, $P(a\leq X\lt b)=P(a\lt X\leq b)=P(a\leq X \leq b)=P(a\lt X\lt b)$, because the probability of $X=c$ is exactly zero; the probability is spread out over an infinite number of values $c$.



- **Cumulative Density Function (CDF)**: A function giving the probability that a continuous random variable is less than or equal to a specific value, $F_X(x) = P(X \leq x)$.

The CDF tells you probabilities on "half intervals" $(\infty, x]$. But it turns out that if you know probabilities on "half intervals," you know everything about the distribution itself. First, notice we automatically have probabilities on any interval $(a, b)$:

```math
P(a\lt X\lt b)=\int_a^bf_X(x) dx=\int_\infty^af_X(x) dx-\int_\infty^bf_X(x) dx=F_X(a)-F_X(b)
```

**Observation:** We have $F_X(x)- F_X(c)=\int_c^xf_X(t) dt$. The Fundamental Theorem of Calculus says that $\frac{d}{dx}\left[\int_c^xf_X(t) dt\right]=f_X(x)$. So,

```math
\frac{d}{dx}\left[F_X(x)- F_X(c)\right]=\frac{d}{dx}\left[\int_c^xf_X(t) dt\right]=f_X(x).
```

The RHS is independent of $c$, so the LHS is, too. (Let $c\to -\infty$.) This shows:

**Theorem:** The PDF $f_X(t)$ and CDF $F_X(t)$ are related by $\frac{d}{dt}\left[F_X(t)\right]=f_X(t)$.



- **Survival Function**: A function giving the probability that a random variable *exceeds* a specific value, $S_X(x) = P(X > x) = 1 - F(x)$.
- **Hazard Function**: A function describing the instantaneous rate at which an event occurs at time $t$, given that the event has not occurred up to time $t$:

```math
h_X(t):=\lim_{\Delta t\to0}\frac{P(t\leq X \leq t+\Delta t \mid X\geq t)}{\Delta t}.
```

**Theorem:** $h_X(t)=\frac{f_X(t)}{S_X(t)}$.

**Proof:**

$h_X(t):=\lim_{\Delta t\to0}\frac{P(t\leq X \leq t+\Delta t \mid X\geq t)}{\Delta t}$

$\hspace{2cm}=\lim_{\Delta t\to0}\frac{P(t\leq X \leq t+\Delta t)}{P(X\geq t) \Delta t}$

$\hspace{2cm}=\frac{\lim_{\Delta t\to0}\frac{P(t\leq X \leq t+\Delta t)}{\Delta t}}{P(X\geq t)}$

$\hspace{2cm}=\frac{\frac{dF_X}{dt}(t)}{S_X(t)}$

$\hspace{2cm}=\frac{f_X(t)}{S_X(t)}\hspace{5cm}\Box$


# 0.1: Background: Poisson Processes and the Exponential Distribution

We will work up to the proof of the following theorem.

> **Theorem:** The only continuous memoryless probability distribution is the exponential distribution.

## Poisson and Exponential Distributions

First recall some definitions.

- **Exponential Distribution**: A continuous probability distribution describing the time between events in a Poisson process, with PDF $f(x) = \lambda e^{-\lambda x}$ for $x \geq 0$.  If an r.v. $X$ has exponential distribution with rate $\lambda$, we write $X \sim \text{Exp}(\lambda)$.
- **Poisson Process**: A stochastic process counting events over time, with events occurring independently at a constant or time-varying rate $\lambda(t)$.
- **Poisson Distribution**: A *discrete* probability distribution giving the probability of $k$ events occurring in a fixed interval, with PDF $F_X(k) = \frac{\lambda^k e^{-\lambda}}{k!}$.

## Memorylessness

 - The distribution of a positive continuously distributed r.v. $X$ is **memoryless** if $P(X>t) = P(X > s+t \mid X > s)$ for all $s,  t\geq 0$.

In words, the probability that an event occurs after a time $t$ (from the start) is the same as the probability that it occurs after an additional time $t$ from a time $s$ given that it has not occurred by time $s$. The time we have to wait for an event is independent of the time we have already waited.

**Lemma:** The distribution of a positive continuously distributed r.v. $X$ is memoryless iff $P(0\lt  X \lt  t) = P(s\lt  X \lt  t+s \mid X > s)$.

**Proof:**

$\hspace{0.5cm}P(X>t) = P(X > s+t \mid X > s)$

$\to 1-P(X>t) = 1-P(X > s+t \mid X > s)$


$\to P(0\lt  X \lt  t) = P(s\lt  X \lt  t+s \mid X > s) \hspace{5cm}\Box$



## Main Theorem

We are ready to prove the main theorem of this section.

**Theorem:** The only continuous memoryless probability distribution is the exponential distribution.

**Proof:**  Observe that for a continuous memoryless distribution, $h(t)=\lambda>0$, a constant. Also, since $S(t)=1-F(t)$, we have $\frac{dS}{dt} = -\frac{dF}{dt}$. It follows that

```math
\lambda = h(t) = \frac{\frac{dF}{dt}}{S(t)} = \frac{-\frac{dS}{dt}}{S(t)} \leftarrow\rightarrow  -\lambda S(t) = \frac{dS}{dt} \leftarrow\rightarrow S(t)=e^{-\lambda t}+ c.
```

The initial condition $S(0)=P(X\geq 0) = 1$ gives $c=0$. Since $\lambda =h(t) = \frac{f(t)}{S(t)}$, it follows that the PDF function $f(t) = \lambda e^{-\lambda t} \quad (t\geq 0). \hspace{5cm}\Box$

**Alternative Proof:** The memoryless condition says

$P(X>t) = P(X > s+t \mid X> s) = \frac{P(X > s+t)}{P(X > s)}$

$\to P(X > s+t) = P(X>t)\cdot P(X>s)$

$\to P(X > at) = P(X>t)^a$ for all $a\in \mathbb{N}$.

Similarly, by substituting $\displaystyle t=\frac{\hat{t}}{a}$, we can also show that $P(X > \frac{t}{b}) = P(X>t)^\frac{1}{b}$ for all $b\in \mathbb{N}$. Thus, $P(X > qt) = P(X>t)^q$ for all $q \in \mathbb{Q}$. By continuity,  $P(X > rt) = P(X>t)^r$ for all $r\in \mathbb{R}$.

This says that as a function of $r$, $P(X>r)$ is an exponential function of the form $P(X>1)^r$.  Or, swapping $t$ for $r$, $P(X>t) = e^{-\lambda t}$, where $\lambda = - \ln(P(X>1)) = -\ln(S(1))$. Thus, we have shown that $S(t) = e^{-\lambda t}$ for some $\lambda > 0$. The CDF is thus $F(t) = 1 - S(t) = 1 - e^{-\lambda t}$. The PDF is $f(t) = \frac{dF}{dt}(t) = \lambda e^{-\lambda t}.\hspace{5cm}\Box$


# 0.2: Time Scaling

 - A Poisson process is **inhomogeneous** if the rate $\lambda$ is not constant over time.

If $\lambda = \lambda(t)$ is an integrable function and $X ∼ \text{Exp}(\lambda)$, then the PDF is $f_X(t) = \lambda(t) e^{-\int_0^t \lambda(s) ds}$, and the CDF is $F_X(t) = 1 - e^{-\int_0^t \lambda(s) ds}$. (Verify using the Fundamental Theorem of Calculus that $F(t) = \int_0^tf(s) ds$.) Let $\Lambda(t)$ be the cumulative rate function defined by

$\displaystyle\hspace{2cm}  \Lambda(t) := \int_0^t \lambda(s)  ds$

Then we can write

$\displaystyle\hspace{2cm} f_X(t) = \lambda(t) e^{-\Lambda(t)}$, and

$\displaystyle\hspace{2cm} F_X(t) = 1 - e^{-\Lambda(t)}$.



**Theorem:** Let $\lambda = \lambda(t)$ be some rate function. If $X \sim \text{Exp}(1)$, then $d(X) \sim \text{Exp}(\lambda)$, where $d = \Lambda^{-1}$, the function inverse of $\Lambda$.

**Example:** Suppose $\lambda(t) = 0.7$. Then $\Lambda(t) := \int_0^t 0.7  ds = 0.7t$, and $d(t) = \frac{t}{0.7}$. By the theorem, if $X \sim \text{Exp}(1)$, then $\frac{1}{0.7}X \sim \text{Exp}(0.7)$.

This says that for a Poisson process, if the rate ($\lambda=0.7$) is smaller (than $1$ in this case), the wait time until the next event ($\frac{1}{0.7}X$) is larger.

**Proof:** We perform a change of variables $Y=d(X)$. By definition,

$F_Y(y) := P(Y \leq y) = P(d(X) \leq y)  \overset{\star}{=} P(X \leq \Lambda(y)) =: F_X(\Lambda(y))$.

Now we use the Chain Rule to differentiate:

$f_Y(y) := \frac{d}{dy}\left[ F_Y(y) \right] = \frac{d}{dy}\left[ F_X(\Lambda(y)) \right] = f_X(\Lambda(y)) \cdot \frac{d}{dy}\left[ \Lambda(y) \right]$.

By the Fundamental Theorem of Calculus, $\frac{d}{dy}\left[ \Lambda(y) \right] = \lambda(y)$. Thus, $f_Y(y) = \lambda(y)e^{-\Lambda(y)}$, the PDF of the distribution $\text{Exp}(\lambda)$. This shows that $Y = d(X) \sim \text{Exp}(\lambda)$. $\hspace{3cm}\Box$

The equality marked with a '$\star$' needs justification, because it's not true in general. If $a \lt  b$, what conditions on a function $f$ guarantee that $f(a) \leq f(b)$? If ($\star$) is to hold, $\Lambda(t)$ must satisfy those conditions.

**Answer:** $f(a) \leq f(b)$ clearly holds if $f$ is monotonically increasing (by definition). Since $\lambda(t)>0$, $\Lambda(t) := \int_0^t \lambda(s)  ds$ is in fact positive and monotonically increasing.
