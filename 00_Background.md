# 0.0: Basic Definitions from College Probability and Statistics 101


- **Probability Density Function (PDF)**: A function describing the relative likelihood of a continuous random variable taking a specific value.

For continuous distributions, this definition isn't good enough, because $`P(X=c)=0`$ for all numbers $`c`$! Instead the PDF describes the probability that an event will occur in some arbitrary interval. If $`f_X(x)`$ is the PDF for the distribution of an r.v. $`X`$, then

```math
P(a \lt X\lt b)=\int_a^bf_X(x) dx
```

**Observation:** *If you know everything about a probability distribution on intervals, then you know the distribution itself.*

**Observation:** For continuous probability distributions, $`P(a\leq X\lt b)=P(a\lt X\leq b)=P(a\leq X \leq b)=P(a\lt X\lt b)`$, because the probability of $`X=c`$ is exactly zero; the probability is spread out over an infinite number of values $`c`$.



- **Cumulative Density Function (CDF)**: A function giving the probability that a continuous random variable is less than or equal to a specific value, $`F_X(x) := P(X \leq x) = \int_{-\infty}^x f_X(w) dw`$.

The CDF tells you probabilities on "half intervals" $`(\infty, x]`$. But it turns out that if you know probabilities on "half intervals," you know everything about the distribution itself. First, notice we automatically have probabilities on any interval $`(a, b)`$:

```math
P(a\lt X\lt b)=\int_a^bf_X(x) dx= \int_{-\infty}^{b}f_X(x) dx -  \int_{-\infty}^{a}f_X(x) dx =F_X(b) - F_X(a)
```

**Observation:** We have $`F_X(x)- F_X(c)=\int_c^xf_X(t) dt`$. The Fundamental Theorem of Calculus says that $`\frac{d}{dx}\left[\int_c^xf_X(t) dt\right]=f_X(x)`$. So,

```math
\frac{d}{dx}\left[F_X(x)- F_X(c)\right]=\frac{d}{dx}\left[\int_c^xf_X(t) dt\right]=f_X(x).
```

The RHS is independent of $`c`$, so the LHS is, too. (Let $`c\to -\infty`$.) This shows:

**Theorem:** The PDF $`f_X(t)`$ and CDF $`F_X(t)`$ are related by $`\frac{d}{dt}\left[F_X(t)\right]=f_X(t)`$.



- **Survival Function**: A function giving the probability that a random variable *exceeds* a specific value, $`S_X(x) = P(X > x) = 1 - F(x)`$.
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


With $S_X(x) = P(X > x) = 1 -  P(X \leq x) = F_X(x)$, then we also have 
```math
P(a\lt X\lt b)= F_X(b) - F_X(a) = S_X(a) - S_X(b)
```


# 0.1: Background: Poisson Processes and the Exponential Distribution

We will work up to the proof of the following theorem.

> **Theorem:** The only continuous memoryless probability distribution is the exponential distribution.

## Poisson and Exponential Distributions

First recall some definitions.

- **Exponential Distribution**: A continuous probability distribution describing the time between events in a Poisson process, with PDF $`f(x) = \lambda e^{-\lambda x}`$ for $`x \geq 0`$.  If an r.v. $`X`$ has exponential distribution with rate $`\lambda`$, we write $`X \sim \text{Exp}(\lambda)`$.
- **Poisson Process**: A stochastic process counting events over time, with events occurring independently at a constant or time-varying rate $`\lambda(t)`$.
- **Poisson Distribution**: A *discrete* probability distribution giving the probability of $`k`$ events occurring in a fixed interval, with PDF $`F_X(k) = \frac{\lambda^k e^{-\lambda}}{k!}`$.

## Memorylessness

 - The distribution of a positive continuously distributed r.v. $`X`$ is **memoryless** if $`P(X>t) = P(X > s+t \mid X > s)`$ for all $`s,  t\geq 0`$.

In words, if we start waiting at time $`0`$ for an event to happen, the probability that the time $`X`$ of the first occurence is at least $`t`$ is $`P(X>t)`$.  To say the distribution is **memoryless** means that if the event has not happened by time $`s`$, then the probability that it does not happen in the next $`t`$ units of time (that is, by time $`t+s`$) is still $`P(X>t)`$.  That is, knowing how long we have waited for an event that has not yet occured gives no information aout how much longer we must wait before it occurs.  

An example of a memoryless discrete distribution comes from rolling a fair die.  If it has not come up $`1`$ after $`n`$ rolls, then our prediction for how many additional rolls are needed until the first time it comes up $`1`$ is the same as if we were about to do our first roll.

The so-called gambler's fallacy that if a result has not happened recently its probability is increasing is based on denying the memoryless property of these distributions.

**Lemma:** The distribution of a positive continuously distributed r.v. $`X`$ is memoryless iff $`P(0\lt  X \lt  t) = P(s\lt  X \lt  t+s \mid X > s)`$.

**Proof:**

$`\hspace{0.5cm}P(X>t) = P(X > s+t \mid X > s)`$

$`\to 1-P(X>t) = 1-P(X > s+t \mid X > s)`$


$`\to P(0\lt  X \lt  t) = P(s\lt  X \lt  t+s \mid X > s) \hspace{5cm}\Box`$



## Main Theorem

We are ready to prove the main theorem of this section.

**Theorem:** The only continuous memoryless probability distribution is the exponential distribution.

**Proof:**  Observe that for a continuous memoryless distribution, $`h(t)=\lambda>0`$, a constant. Also, since $`S(t)=1-F(t)`$, we have $`\frac{dS}{dt} = -\frac{dF}{dt}`$. It follows that

```math
\lambda = h(t) = \frac{\frac{dF}{dt}}{S(t)} = \frac{-\frac{dS}{dt}}{S(t)} \leftarrow\rightarrow  -\lambda S(t) = \frac{dS}{dt} \leftarrow\rightarrow S(t)=e^{-\lambda t}+ c.
```

The initial condition $`S(0)=P(X\geq 0) = 1`$ gives $`c=0`$. Since $`\lambda =h(t) = \frac{f(t)}{S(t)}`$, it follows that the PDF function $`f(t) = \lambda e^{-\lambda t} \quad (t\geq 0). \hspace{5cm}\Box`$

**Alternative Proof:** The memoryless condition says

$`P(X>t) = P(X > s+t \mid X> s) = \frac{P(X > s+t)}{P(X > s)}`$

$`\to P(X > s+t) = P(X>t)\cdot P(X>s)`$

$`\to P(X > at) = P(X>t)^a`$ for all $`a\in \mathbb{N}`$.

Similarly, by substituting $`\displaystyle t=\frac{\hat{t}}{a}`$, we can also show that $`P(X > \frac{t}{b}) = P(X>t)^\frac{1}{b}`$ for all $`b\in \mathbb{N}`$. Thus, $`P(X > qt) = P(X>t)^q`$ for all $`q \in \mathbb{Q}`$. By continuity,  $`P(X > rt) = P(X>t)^r`$ for all $`r\in \mathbb{R}`$.

This says that as a function of $`r`$, $`P(X>r)`$ is an exponential function of the form $`P(X>1)^r`$.  Or, swapping $`t`$ for $`r`$, $`P(X>t) = e^{-\lambda t}`$, where $`\lambda = - \ln(P(X>1)) = -\ln(S(1))`$. Thus, we have shown that $`S(t) = e^{-\lambda t}`$ for some $`\lambda > 0`$. The CDF is thus $`F(t) = 1 - S(t) = 1 - e^{-\lambda t}`$. The PDF is $`f(t) = \frac{dF}{dt}(t) = \lambda e^{-\lambda t}.\hspace{5cm}\Box`$


**Another Alternate Proof**
Outline: The intuition behind this proof is based on dividing $`[0,t]`$ into many small intervals of length $`a=t/n`$ for some large integer $n$ and noting that the event does not happen in $`[0,t]`$ iff it does not happen in any subinterval.  We can estimating the probability that an event does not happen in each interval and use this to estimate the probability that the event does not happen in $[0,t]$.  We then take $`a \to 0`$ (or equivalently $`n \to \infty`$) to find $`P(X>t)`$.


Steps:

- Following steps in the previous proof, $`P(X>an) = P(X>a)^n`$ for all $`n \in \mathbb{N}`$. 

- Additionally, if $`P(X>0)=1`$, then $`P(X>a) = 1-\int_0^a f_X(w) dw \approx 1- f_X(0) a`$ where the error in the approximation is $\mathcal{o}(a^2)$

Then for any $`t>0`$, we can write $`t=an`$ where $`n \in \mathbb{N}`$ can be arbitrarily large, and $`a=t/n`$ gets small as $`n`$ grows.

$` P(X>t) = (1 - a f_X(0) + \mathcal{o}(a^2))^n = ([1-a f_X(0)] + \mathcal{o}(a^2))^n`$

$`P(X>t) = (1-af_X(0))^n + \cdots`$

The neglected terms (the $`\cdots`$) are of size $na^2$, $n^2a^4$, etc.  Since $`a = t/n`$, the value of $na^2$ is $t^2/n$.  

We hold $`t`$ fixed, but take the number of subintervals $n \to \infty$. The neglected terms will go to zero and $`(1-\frac{t}{n} f_X(0))^n \to e^{-t f_X(0)}`$.

Thus 

$`P(X>t) = e^{-\lambda t}`$ where $`\lambda = f_X(0)`$.$`\hspace{5cm}\Box`$


# 0.2: Time Scaling

 #### A Poisson process is **inhomogeneous** if the rate $`\lambda`$ is not constant over time.

We first present an explanation of the intuition underlying the idea of time scaling of inhomogeneous Poisson processes.

For a **homogeneous** Poisson process, the time to the next event has an exponential distribution whose rate is the rate of the process.  We will be able to rescale our measurement time so that an inhomogeneous Poisson process also has an exponential distribution in our rescaled time.  Intuitively, the way we do this is to consider a clock that ticks at rate $\lambda(t)$, with the clock's measured time passing faster when the rate is higher.  For this clock, the measured "time" to the first event is exponentially distributed.  

This variable clock's measured time is $`\int_0^t \lambda(s) ds`$, the cumulative hazard.  The cumulative hazard when the first event occurs is exponentially distributed with rate $1$.  The time of that event is the time at which the cumulative hazard reaches the threshold.

Let's reframe this in the context of an infectious disease.  The time until an individual gets infected is a random variable.  If we measure in terms of the cumulative exposure to the pathogen, the cumulative exposure $`Y`$ measured by the variable clock at the first successful transmission is exponentially distributed with constant rate.  However, since the exposure rate itself varies in time (based on behavior, amount of infection present, etc), to calculate the time of infection, we need to determine the real time $`X`$ when the variable clock reaches $`Y`$, that is the time $`X`$ at which $`\Lambda(X) = \int_0^X \lambda(s) ds`$ equals $`Y`$.

Now let's reframe this mathematically.  We are focused on the distribution of the waiting time $`X`$ until the first event.

If $`\lambda = \lambda(t)`$ is an integrable function and $`X ∼ \text{Exp}(\lambda)`$, then the PDF of $`X`$ is $`f_X(t) = \lambda(t) e^{-\int_0^t \lambda(s) ds}`$, and the CDF is $`F_X(t) = 1 - e^{-\int_0^t \lambda(s) ds}`$. (Verify using the Fundamental Theorem of Calculus that $`F(t) = \int_0^tf(s) ds`$.) Let $`\Lambda(t)`$ be the cumulative rate function defined by

$`\displaystyle\hspace{2cm}  \Lambda(t) := \int_0^t \lambda(s)  ds`$

Then we can write

$`\displaystyle\hspace{2cm} f_X(t) = \lambda(t) e^{-\Lambda(t)}`$, and

$`\displaystyle\hspace{2cm} F_X(t) = 1 - e^{-\Lambda(t)}`$.



**Theorem:** Let $`\lambda = \lambda(t)`$ be some rate function. If $`Y \sim \text{Exp}(1)`$, then $`X=d(Y) \sim \text{Exp}(\lambda)`$, where $`d = \Lambda^{-1}`$, the function inverse of $`\Lambda`$.

In the variable clock analogy, if the measured time $`Y`$ to the first event as measured by the variable clock is exponentially distributed with rate $1$, then the real time of the first event is $`X=\Lambda^{-1}(Y)`$ and it has distribution $`X \sim \text{Exp}(\lambda)`$

**Example:** Suppose $`\lambda = 0.7`$ is homogeneous. Then $`Y=\Lambda(X) := \int_0^X 0.7  ds = 0.7X`$, and $`X=d(Y) = \frac{Y}{0.7}`$. By the theorem, if $`Y \sim \text{Exp}(1)`$, then $`X=\frac{1}{0.7}Y \sim \text{Exp}(0.7)`$.

This says that for a Poisson process, if the rate ($`\lambda=0.7`$) is smaller (than $`1`$ in this case), the wait time until the next event ($`\frac{1}{0.7}X`$) is larger.


**Proof:** We perform a change of variables $`X=d(Y)`$. By definition,

$`F_X(x) := P(X \leq x) = P(d(Y) \leq x)  \overset{\star}{=} P(Y \leq \Lambda(x)) =: F_Y(\Lambda(x))`$.

Now we use the Chain Rule to differentiate:

$`f_X(x) := \frac{d}{dx}\left[ F_X(x) \right] = \frac{d}{dx}\left[ F_Y(\Lambda(x)) \right] = f_Y(\Lambda(x))) \cdot \frac{d}{dx}\left[ \Lambda(x) \right]`$.

By the Fundamental Theorem of Calculus, $`\frac{d}{dx}\left[ \Lambda(x) \right] = \lambda(x)`$. Thus, $`f_X(x) = \lambda(x)e^{-\Lambda(x)}`$, the PDF of the distribution $`\text{Exp}(\lambda)`$. This shows that $`X = d(Y) \sim \text{Exp}(\lambda)`$. $`\hspace{3cm}\Box`$


The equality marked with a '$`\star`$' needs justification, because it's not true in general. If $`a \lt  b`$, what conditions on a function $`f`$ guarantee that $`f(a) \leq f(b)`$? If ($`\star`$) is to hold, $`\Lambda(t)`$ must satisfy those conditions.

**Answer:** $`f(a) \leq f(b)`$ clearly holds if $`f`$ is monotonically increasing (by definition). Since $`\lambda(t)>0`$, $`\Lambda(t) := \int_0^t \lambda(s)  ds`$ is in fact positive and monotonically increasing.


**Another proof:**

Let $\lambda = \lambda(t)$ be some non-negative real-valued rate function of a Poisson point process, i.e. the hazard function. When $\lambda(t)$ is not constant, the Poisson point process is inhomogeneous. 

Let $c(t)$ be the cumulative rate function, s.t.

```math
c(t) = \int_0^t \lambda(s)ds 
```
With $\lambda(t)$ a non-negative real-valued function, the cumulative of this must always be monotonically increasing. 

$\lambda(t)$ as the hazard function also satisfies the relationship,

```math
\lambda(t) = \frac{f(t)}{S(t)} = \frac{1}{S(t)} \frac{dF}{dt} = \frac{-1}{S(t)}\frac{dS}{dt} = -\frac{d}{dt}\Big(lnS(t)\Big)
```

Then 

```math
ln(S(t)) = -\int_0^t \lambda(s)ds = -c(t)
```

Therefore, we can write the survival function $S(t) = P(T > t)$ as

```math
S(t) = e^{-c(t)}
```

Now since $y=c(t)$ is a monotonically increasing function, this means that it is also invertible, i.e., there is another function $d(y)$ s.t. $d(c(t)) = t$ and $c(d(y)) = y$.

Therefore, $c(T) > y \Leftrightarrow T > d(y)$ since c(t) and d(y) uniquely map to each other. 

Then

```math
P(c(T) > y) = P(T > d(y)) = S(d(y)) = e^{-c(d(y))} = e^{-y} = S^{*}(y)
```

where $S^{*}(y)$ is the survival function of events in the cumulative space. 

This means that $c(T) = y \sim Exp(1)$, i.e. inter event distances in the cumulative space can be sampled with an exponential distribution of rate equal to $1$.  $`\hspace{3cm}\Box`$


Now, what does this mean for us in terms of being able to sample the time to the next infection when the rate process is inhomogeneous?

Since we can sample the inter event distances in the cumulative space with $\sim Exp(1)$, then we can get the cumulative values $y_i$ by summing the sampled inter event values $\Delta y_i$, and invert each $y_i$ with $d(y)$ to get the time for each event $t_i$. 

Algorithmically, this means if we want to identify the times of a sequence of events from an inhomogeneous Poisson process, we can 

1. Sample the scaled inter-interval "times" $\Delta y_i$ from $Exp(1)$
2. Calculate the scaled times $y_i$ as $y_i = \Sigma_{i' = 1}^i \Delta y_{i'} $
3. Invert $y_i$ to get the true times $t_i$ through $t_i = d(y_i)$
4. Calculate the inter event times $\Delta t_i = t_i - t_{i-1}$

where $\Delta t_1 = t_1 - t_0 = t_1 - 0 = t_1$ since the time of the zeroeth event would be 0.



