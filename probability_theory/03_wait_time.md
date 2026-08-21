# Modeling Wait Time

This chapter studies a random variable that models a wait time, that is, time till some event.

## Specialization to Wait Times

In the previous chapter we assembled a cast of characters — the CDF $`F_X`$, the survival function $`S_X`$, the density $`f_X`$, the hazard $`h_X`$, and the cumulative hazard $`H_X`$ — and worked out how they are related. That development was deliberately distribution-agnostic: nothing there assumed $`X`$ meant anything in particular. We now specialize to the case that $`X \ge 0`$ is a random variable representing a *wait time*: the time until some event occurs. The event might be the failure of a component, the death of an organism, the decay of an atom, or the arrival of the next bus. Because a wait time is never negative, we take $`f_X`$ and $`F_X`$ to vanish for $`x < 0`$, and we adopt the convention
```math
H_X(x) = \int_0^x h_X(w)\,dw \qquad (x \ge 0)
```
for the lower limit throughout, dropping the $`-\infty`$ that appeared in the general theory.


Recall each object and its wait-time interpretation.

- The *CDF* $`F_X(x) = P(X \le x)`$ is the probability that the event has already happened by time $`x`$. It climbs from $`0`$ toward $`1`$ as we give the event more and more time to occur.

- The *survival function* $`S_X(x) = 1 - F_X(x) = P(X > x)`$ is the probability that we are still waiting at time $`x`$ — the event has not yet happened. It is the mirror image of $`F_X`$, starting at $`1`$ and decaying toward $`0`$.

- The *density* $`f_X = F_X'`$ is the unconditional instantaneous event-rate: recalling
  ```math
  f_X(x) = \lim_{h \to 0^+} \frac{P(x < X \le x+h)}{h},
  ```
  it measures how much probability is concentrated near time $`x`$, as weighed against the entire original population of possible outcomes. It gives the fraction of everyone who started waiting whose event occurs right around time $`x`$.

- The *hazard* $`h_X`$ will be a primary object of our study this chapter, so we discuss it in detail below.

- The *cumulative hazard* $`H_X(x) = \int_0^x h_X(w)\,dw`$ is the total risk accumulated by time $`x`$, the running tally of exposure a survivor has absorbed just by lasting that long.

### The hazard

The density $`f_X`$ measures risk against everyone, but that is usually not the question we need to answer. If I am holding a component that has already run without failure for a year, I do not care about the failure rate averaged over all components including the ones that died in week one. I care about my component's risk, given that it is still running. That conditional quantity is the hazard.

**Def (recalled).** The *hazard* (hazard rate, hazard function) of $`X`$ is
```math
\begin{aligned}
h_X(x) &= \lim_{h \to 0^+} \frac{1}{h}\, P\big(x < X \le x+h \;\big|\; X > x\big) \\[4pt]
       &= \lim_{h \to 0^+} \frac{1}{h}\cdot\frac{P(x < X \le x+h)}{P(X > x)} \\[4pt]
       &= \frac{f_X(x)}{S_X(x)}.
\end{aligned}
```

The density gives the instantaneous event rate relative to the whole population; the hazard gives the instantaneous event rate relative only to those still waiting at $`x`$ — the survivors. Dividing by $`S_X(x)`$ is exactly the act of throwing away everyone whose event already happened and re-normalizing among the ones who are left. So

> $`h_X(x)`$ is the instantaneous rate of the event at time $`x`$ for an individual who has made it to time $`x`$.


### The hazard determines the distribution

In Chapter 1 we derived the relationships among these objects and, in particular, the identity
```math
S_X(x) = e^{-H_X(x)} \quad\longleftrightarrow\quad H_X(x) = -\log S_X(x).
```

This identity lets us reverse the usual order of calculation. In a first course one is handed a density $`f_X`$, integrates to get $`F_X`$, subtracts to get $`S_X`$, and only then divides to recover the hazard. The density is treated as primary and the hazard as a derived curiosity.

For wait times we will run the machine in the opposite direction. The hazard is the thing we can reason about mechanistically — it is a statement about the physics of the process, the instantaneous risk faced by a survivor. Given a hazard, the identity above hands us the entire distribution for free:
```math
h_X \;\xrightarrow{\;\int_0^x\;}\; H_X \;\xrightarrow{\;e^{-(\cdot)}\;}\; S_X \;\xrightarrow{\;1 - (\cdot)\;}\; F_X \;\xrightarrow{\;\frac{d}{dx}\;}\; f_X.
```
Every object downstream of $`h_X`$ is determined by it. Thus, instead of choosing a wait-time distribution directly, we can specify how the conditional event rate changes with elapsed time.

## The Homogeneous Case: Constant Hazard

The simplest wait-time model has a hazard that is constant in time. Under this assumption, the instantaneous event rate for an individual still waiting does not depend on how long that individual has waited.

**Def:** A wait-time model is *homogeneous* if its hazard is constant in time:
```math
h_X(x) \equiv \lambda
```
for some fixed $`\lambda > 0`$.

The name records the fact that the process looks the same at every moment — homogeneous in time. The value of $`h_X(x)`$ does not depend on $`x`$: an individual who has waited a long time faces the same instantaneous event rate as one who has just begun waiting.

### Deriving the exponential distribution

Substituting the constant hazard into the cumulative-hazard formula gives
```math
H_X(x) = \int_0^x h_X(w)\,dw = \int_0^x \lambda\,dw = \lambda x,
```
and the identity $`S_X = e^{-H_X}`$ delivers the survival function immediately:
```math
S_X(x) = e^{-\lambda x}.
```
Everything else follows by the relationships we already have. The CDF is the complement of survival,
```math
F_X(x) = 1 - S_X(x) = 1 - e^{-\lambda x},
```
and the density is its derivative,
```math
f_X(x) = F_X'(x) = \lambda e^{-\lambda x},
```
each understood to hold for $`x \ge 0`$ and to be zero for $`x < 0`$.

**Def:** A wait time with density $`f_X(x) = \lambda e^{-\lambda x}`$ on $`x \ge 0`$ is *exponentially distributed* with rate $`\lambda`$, written $`X \sim \text{Exponential}(\lambda)`$.

Thus, a constant hazard $`h_X(x) = \lambda`$ determines the exponential distribution with rate $`\lambda`$.

### The rate parameter and mean wait time

The constant $`\lambda`$ is the hazard itself, so it inherits the hazard's interpretation directly: it is the instantaneous rate of the event among those still waiting, and here that rate is the same at every elapsed time. Larger $`\lambda`$ means a survivor is in more danger at every instant, so the event tends to come sooner and the survival curve $`e^{-\lambda x}`$ falls away more steeply.

The parameter $`\lambda`$ has units of inverse time, and its reciprocal gives the mean wait time. For $`X \sim \text{Exponential}(\lambda)`$,
```math
E[X] = \int_0^\infty x\,\lambda e^{-\lambda x}\,dx = \frac{1}{\lambda},
```
which matches the interpretation: an event occurring at conditional rate $`\lambda`$ per unit time has mean wait time $`1/\lambda`$. A component with hazard $`\lambda = 0.01`$ per hour has a mean lifetime of $`100`$ hours.

### A preview of memorylessness

Constant hazard also determines the distribution of the remaining wait time. In particular, elapsed waiting time does not change that distribution. This memoryless property characterizes the exponential distribution among continuous wait-time distributions; we return to it after considering time-varying hazards.

## The Inhomogeneous Case: Time-Varying Hazard

A constant hazard is the exception rather than the rule. Usually a survivor's risk depends on how long it has already waited: a machine part wears and grows more likely to fail, while a new device may survive an initial fragile period and settle into safety.

The pipeline from the first section already handles this. Nothing in it required the hazard to be constant; we chose a constant hazard only to see the simplest case first. We now allow the hazard to be any nonnegative function of time.

**Def:** A wait time $`X`$ is *inhomogeneous* if its hazard varies in time:
```math
h_X(x) = \lambda(x)
```
for some function $`\lambda(x) \ge 0`$.

We continue to write the hazard with the letter $`\lambda`$, now carrying an argument, so that the constant-hazard case sits inside this one as the special choice $`\lambda(x) \equiv \lambda`$. Here $`\lambda(x)`$ is the hazard: the instantaneous risk of the event at time $`x`$ for an individual still waiting at time $`x`$.

The cumulative hazard is again the integral of the hazard, though the integral no longer collapses to a product:
```math
H_X(x) = \int_0^x h_X(w)\,dw = \int_0^x \lambda(w)\,dw.
```
The identity $`S_X = e^{-H_X}`$ gives the survival function,
```math
S_X(x) = \exp\!\left(-\int_0^x \lambda(w)\,dw\right),
```
and the CDF and density follow as before:
```math
\begin{align*}
F_X(x) &= 1 - \exp\!\left(-\int_0^x \lambda(w)\,dw\right), \\
f_X(x) &= \lambda(x)\,\exp\!\left(-\int_0^x \lambda(w)\,dw\right).
\end{align*}
```

The density factors into two pieces with direct readings. The exponential term is $`S_X(x)`$, the probability of surviving to $`x`$. The prefactor $`\lambda(x)`$ is the hazard at $`x`$, the risk faced right then by a survivor. Their product is the unconditional density: to have the event near $`x`$, an individual must survive to $`x`$ and then succumb, and $`f_X(x) = \lambda(x)\,S_X(x)`$ records those two requirements. This is the relationship $`h_X = f_X/S_X`$ from the first section, solved for $`f_X`$.

Setting $`\lambda(x) \equiv \lambda`$ recovers $`\int_0^x \lambda\,dw = \lambda x`$ and returns the exponential of the previous section. The exponential is the member of this family whose hazard is flat.

**Note:** For $`S_X`$ to describe an event that is certain to occur eventually, the total accumulated hazard must diverge:
```math
\int_0^\infty \lambda(w)\,dw = \infty.
```
If the integral stayed finite, $`S_X(x)`$ would level off at a positive value as $`x \to \infty`$, leaving a nonzero probability that the event never happens.

**Note:** The hazard $`\lambda(x)`$ is a rate, not a probability. It carries units of inverse time and may exceed $`1`$; only $`F_X`$ and $`S_X`$ are confined to $`[0,1]`$.

With the inhomogeneous case in hand, we can build a wait-time distribution from any hazard we can justify, and the exponential sits at the single point where that hazard is constant. The next section takes up the property that distinguishes this point from all the others.

## Memorylessness

A constant hazard means a survivor's instantaneous risk does not depend on how long it has waited. This abstract property is called *memorylessness*, because a *memoryless* distribution does not "remember" any previous time spent waiting. We will characterize all continuous distributions having this property and discover the surprising conclusion: the exponential distribution is the only one!

Suppose a component has run without failure up to time $`s`$, and we ask how much longer it will last. The remaining wait is $`X - s`$, and we want its distribution given that the component has survived to $`s`$, namely $`P(X - s > t \mid X > s)`$. A wait time is *memoryless* when this conditional distribution of remaining time does not depend on $`s`$: a survivor's future looks the same no matter how long it has already waited.

**Def:** A wait time $`X`$ is *memoryless* if
```math
P(X > s + t \mid X > s) = P(X > t)
```
for all $`s, t \ge 0`$.

The right side is the survival function of a fresh arrival evaluated at $`t`$. The left side is the probability that a survivor at time $`s`$ lasts at least $`t`$ longer. Their equality says that a used component in working order is indistinguishable from a new one.

### Memorylessness in terms of the survival function

We can rewrite this defintion in terms of the survival function. Using $`P(X > s+t \mid X > s) = P(X > s+t)\,/\,P(X > s)`$ and $`S_X(x) = P(X > x)`$, memorylessness becomes
```math
\frac{S_X(s+t)}{S_X(s)} = S_X(t),
```
or, clearing the denominator,
```math
S_X(s+t) = S_X(s)\,S_X(t) \qquad \text{for all } s, t \ge 0.
```
The probability of surviving $`s + t`$ factors into the probability of surviving the first $`s`$ and the probability of surviving a further $`t`$, with the second factor unaffected by the first. This is the functional form we work with.

### Characterizing memoryless distributions

Memorylessness is now a question about a single functional equation: which continuous survival functions satisfy $`S_X(s+t) = S_X(s)\,S_X(t)`$? The exponential is one solution, since $`e^{-\lambda(s+t)} = e^{-\lambda s}\,e^{-\lambda t}`$ by the law of exponents. The following theorem shows it is the only one.

**Theorem.** Let $`X`$ be a continuous wait time with $`P(X > 0) = 1`$. If $`X`$ is memoryless, then $`S_X(x) = e^{-\lambda x}`$ for some $`\lambda > 0`$; that is, $`X`$ is exponentially distributed.

**Proof.** Write $`g(x) := S_X(x)`$. Memorylessness gives the functional equation
```math
g(s + t) = g(s)\,g(t) \qquad (s, t \ge 0).
```
We solve this equation first in general, using only that $`g`$ is continuous and not identically zero, and afterward specialize to the survival function of a wait time.

*Solving the functional equation.* Suppose $`g : [0,\infty) \to \mathbb{R}`$ is continuous, not identically zero, and satisfies $`g(s+t) = g(s)\,g(t)`$ for all $`s, t \ge 0`$. We claim $`g(x) = c^x`$ for a constant $`c > 0`$.

Setting $`s = t = 0`$ gives $`g(0) = g(0)^2`$, so $`g(0) \in \{0, 1\}`$. If $`g(0) = 0`$, then $`g(x) = g(x)\,g(0) = 0`$ for every $`x`$, making $`g`$ identically zero; since it is not, $`g(0) = 1`$.

Let $`c := g(1)`$. Applying the functional equation repeatedly, for any positive integer $`n`$,
```math
g(n) = g(\underbrace{1 + \cdots + 1}_{n}) = g(1)^n = c^n.
```
The same additive splitting applied to $`1 = \tfrac1n + \cdots + \tfrac1n`$ gives $`g(1) = g(\tfrac1n)^n`$, so $`g(\tfrac1n) = c^{1/n}`$. Combining these, for any positive integers $`m, n`$,
```math
g\!\left(\tfrac{m}{n}\right) = g\!\left(\tfrac1n\right)^m = c^{m/n},
```
so $`g(q) = c^q`$ for every nonnegative rational $`q`$. Every real $`x \ge 0`$ is a limit of rationals $`q_k \to x`$, and by continuity
```math
g(x) = \lim_{k \to \infty} g(q_k) = \lim_{k \to \infty} c^{q_k} = c^x.
```
Finally, $`c > 0`$: if $`c = g(1)`$ were zero, then $`g(\tfrac1n)^n = 0`$ would force $`g(\tfrac1n) = 0`$, and continuity at $`0`$ would give $`g(0) = 0`$, contradicting $`g(0) = 1`$. This proves $`g(x) = c^x`$ with $`c > 0`$.

*Specializing to a survival function.* The survival function $`S_X`$ is continuous by hypothesis and not identically zero, since $`S_X(0) = 1`$, so $`S_X(x) = c^x`$ for some $`c > 0`$. Because $`X`$ eventually occurs, $`S_X`$ is not identically $`1`$, forcing $`c < 1`$. With $`0 < c < 1`$ we may write $`c = e^{-\lambda}`$ for a unique $`\lambda > 0`$, giving
```math
S_X(x) = c^x = e^{-\lambda x}.
```
This is the survival function of $`\text{Exponential}(\lambda)`$. $`\quad\blacksquare`$

**Note:** The functional equation $`g(s+t) = g(s)\,g(t)`$ is one of Cauchy's functional equations. Without a regularity assumption such as continuity it admits pathological solutions, built with a Hamel basis for $`\mathbb{R}`$ over $`\mathbb{Q}`$, that are not of the form $`c^x`$. Continuity of $`S_X`$ rules these out, which is why the theorem requires a continuous wait time. The monotonicity of a survival function would serve equally well in place of continuity.

#### An easier proof in the differentiable case

If we assume in addition that $`S_X`$ is differentiable, there is a much shorter proof. Start with the functional equation $`S(s+t) = S(s)\,S(t)`$ and differentiate both sides with respect to $`t`$:

```math
S'(s+t) = S(s)\,S'(t).
```

Evaluating at $`t = 0`$,

```math
S'(s) = S'(0)\,S(s).
```

Since $`S`$ is nonincreasing, $`S'(0) \le 0`$; write $`-\lambda := S'(0)`$. This is the first-order linear equation

```math
S'(s) = -\lambda\,S(s), \qquad S(0) = 1,
```

whose unique solution is $`S(s) = e^{-\lambda s}`$. To justify our claim of uniqueness, we could cite the standard existence-uniqueness theorem from ordinary differential equations, or we can simply observe that from $`S' = -\lambda S`$ we have $`\frac{d}{ds}\big(e^{\lambda s} S(s)\big) = e^{\lambda s}\big(S'(s) + \lambda S(s)\big) = 0`$, so $`e^{\lambda s} S(s)`$ is constant, equal to its value $`1`$ at $`s = 0`$.

Of course, our goal was to characterize all memoryless continuous distributions, not just the differentiable ones, so we needed the more powerful proof.
