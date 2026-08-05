# From Wait Time to Poisson Process

In this chapter we will study the concept of a Poisson process from two directions at once: we will construct a Poisson process from the bottom up by stacking wait times; and we will show that our bottom-up construction is precisely the object characterized by the top-down axiomatic definition of Poisson process.

## Stacking Wait Times

So far a single random variable $`X`$ has carried one event: the wait until it happens, once. Of course many processes do not stop after one event: a detector registers one particle, then waits for the next; a server handles a request, then waits for another, and so on. To describe a stream of events rather than a single one, we stack wait times end to end. Let $`\{X_n\}_{n=1}^\infty`$ be a sequence of wait times, where $`X_n`$ is the time between the $`(n-1)`$-th event and the $`n`$-th. We take the $`X_n`$ to be independent and identically distributed, each with the same distribution as a single wait time $`X \ge 0`$. The first event happens after waiting $`X_1`$; the second happens $`X_2`$ later; and so on.

**Def:** The *arrival times* $`\{T_n\}_{n=0}^\infty`$ are the partial sums of the wait times,
```math
T_0 := 0, \qquad T_n := X_1 + X_2 + \cdots + X_n = \sum_{k=1}^n X_k,
```
so that $`T_n`$ is the time at which the $`n`$-th event occurs.

Each $`T_n`$ marks a point on the time axis where an event lands. Because every $X_k$ is nonnegative, the arrival times are nondecreasing, $0 = T_0 \le T_1 \le T_2 \le \cdots$. If in addition the wait times are continuous, then each $X_k > 0$ with probability one, and the arrival times are almost surely strictly increasing. The sequence $`\{T_n\}_{n=0}^\infty`$ is one realization of the event stream: a list of the instants at which something happened.

## The Counting Process

The arrival times record *when* events happen. The complementary view records *how many* have happened by a given time. This is the object we ultimately want, because it turns a stream of events into a function of time we can evaluate at any instant.

**Def:** The *counting process* $`\{N(t)\}_{t \ge 0}`$ associated with the arrival times $`\{T_n\}_{n=0}^\infty`$ is
```math
N(t) := \max\{\, n : T_n \le t \,\},
```
the number of events that have occurred in the interval $`(0, t]`$.

For each fixed $`t`$, the count $`N(t)`$ is a random variable — it depends on the random wait times through the arrival times. Viewed as $`t`$ varies, $`\{N(t)\}_{t \ge 0}`$ is a family of random variables indexed by time, a *stochastic process*. A single realization of the wait times produces a single sample path: a function of $`t`$ that starts at $`N(0) = 0`$ and jumps up by one at each arrival time, holding constant in between. The path is a right-continuous staircase, flat on each interval $`[T_n, T_{n+1})`$ and stepping up by one at $`T_{n+1}`$.

We now have three descriptions of the same event stream, each carrying the same information and each suited to a different question:

- the wait times $`\{X_n\}_{n=1}^\infty`$ — how long between consecutive events;
- the arrival times $`\{T_n\}_{n=0}^\infty`$ — when each event occurs;
- the counting process $`\{N(t)\}_{t \ge 0}`$ — how many events have occurred by time $`t`$.

Each determines the others. The arrival times are partial sums of the wait times, the wait times are successive differences of the arrival times, and the counting process is read off the arrival times by the definition above.

Translating between counts and arrival times is especially useful for computations. The basic translation is

```math
N(t) \ge n \quad\Longleftrightarrow\quad T_n \le t \qquad (n \ge 1,\ t \ge 0),
```

which holds because "at least $`n`$ events have occurred by time $`t`$" and "the $`n`$-th event has occurred by time $`t`$" describe the same situation.

In particular, taking probabilities of both sides gives

```math
P\big(N(t) \ge n\big) = P\big(T_n \le t\big).
```

The left-hand side cannot be addressed directly, whereas we can find the distribution of the right-hand side, since it is the sum of independent wait times.

## The Poisson Process, Defined

The previous section built a concrete object from the ground up: a counting process $`\{N(t)\}_{t \ge 0}`$ assembled from a sequence of independent, identically distributed wait times. If we take those wait times to be exponential — the memoryless wait time singled out in the first half of this chapter — the process we have constructed is exactly the object characterized by the standard mathematical definition of a Poisson process. We state that definition now, in full, and then use it as the roadmap for the rest of the chapter: each clause is a property we will prove our construction satisfies.

A Poisson process is a particular kind of stochastic process.

**Def:** A *stochastic process* is an indexed family of random variables $`\{X_t\}_{t \in T}`$, all defined on a common probability space $`(\Omega, \mathcal{F}, P)`$, where $`T`$ is some index set.

For a *time-indexed* process, $`T = [0, \infty)`$ (or sometimes $`\mathbb{R}`$), and $`X_t`$ represents the state of some system at time $`t`$.

A stochastic process can be viewed two ways:

- **as a family of r.v.s**: for fixed $`t`$, $`X_t : \Omega \to \mathbb{R}`$ is an ordinary random variable.
- **as a random function**: for fixed $`\omega \in \Omega`$, the map $`t \mapsto X_t(\omega)`$ is a deterministic function of $`t`$, called a *sample path* (or *realization*) of the process.

The content of a stochastic process lies not just in the marginal distribution of each $`X_t`$, but in the joint structure relating $`X_s`$ and $`X_t`$ across different times.

The formal definition of a Poisson process makes use of the following discrete probability distribution, although the reasons will not be immediately clear.

**Def:** The *Poisson distribution with rate $`\lambda`$*, denoted $`\text{Poisson}(\lambda)`$, is a discrete distribution defined to have probability mass function $`f_N(k) := \frac{\lambda^{k} e^{-\lambda}}{k!}`$ and cumulative distribution function $`F_N(k) := e^{-\lambda} \sum_{j=0}^{\lfloor k \rfloor} \frac{\lambda^{j}}{j!}`$.

**Def:** Let $`\lambda : [0,\infty) \to [0,\infty)`$ be a (measurable, locally integrable) *rate function*. A *Poisson process* with rate $`\lambda(t)`$ is a stochastic process $`\{N(t)\}_{t \geq 0}`$, where $`N(t)`$ is a random variable representing the number of events that have occurred in $`(0, t]`$, satisfying:

1. $`N(0) = 0`$ — count starts at zero at the origin.
2. $`\{N(t)\}_{t\geq0}`$ has independent increments: for any $`0 \le t_0 < t_1 < \cdots < t_k`$, the random variables $`N(t_1)-N(t_0),\ N(t_2)-N(t_1),\ \ldots,\ N(t_k)-N(t_{k-1})`$ are independent.
3. For $`0 \le s < t`$, $`N(t) - N(s) \sim \text{Poisson}\!\left(\int_s^t \lambda(u)\,du\right)`$.

The definition specifies a Poisson process by its properties, saying nothing about how to build one. Our construction runs the other way: it produces a specific process and leaves its properties to be discovered. The rest of the chapter connects the two ends by proving that the constructed process meets each clause of the definition, and the clauses tell us in what order to proceed.

## The Distribution of the $`n`$-th Arrival

From here on we specialize the construction to exponential wait times. Let $`\{X_n\}_{n=1}^\infty`$ be independent, each $`X_n \sim \text{Exponential}(\lambda)`$ for a fixed rate $`\lambda > 0`$, and let $`\{T_n\}_{n=0}^\infty`$ be the arrival times $`T_n = X_1 + \cdots + X_n`$. This is the homogeneous case, where every gap between events is drawn from the same memoryless distribution.

To reach the distribution of the count $`N(t)`$ we first need the distribution of a single arrival time $`T_n`$. The counting identity of the previous section, $`P(N(t) \ge n) = P(T_n \le t)`$, makes this the natural first step: once we know how $`T_n`$ is distributed, the count follows. And $`T_n`$ is a sum of independent random variables whose common distribution we already understand, so it is within reach.

Recall that for independent $`Y`$ and $`Z`$ with densities $`f_Y`$ and $`f_Z`$, the density of $`Y + Z`$ is the convolution of their densities:
```math
f_{Y+Z}(t) = \int_{-\infty}^\infty f_Y(s)\, f_Z(t - s)\, ds.
```

Each arrival time is the previous one plus an independent wait, $`T_{n+1} = T_n + X_{n+1}`$, so convolution lets us pass from the density of $`T_n`$ to the density of $`T_{n+1}`$. Because all the variables here are nonnegative, the densities vanish for negative arguments, and the integral collapses to the range where both factors are supported: $`f_{T_n}(s)`$ requires $`s \ge 0`$, and $`f_X(t - s)`$ requires $`s \le t`$. The convolution therefore runs only over $`0 \le s \le t`$,
```math
f_{T_{n+1}}(t) = \int_0^t f_{T_n}(s)\, f_X(t - s)\, ds \qquad (t \ge 0),
```
and is zero for $`t < 0`$.

We now compute $`f_{T_n}`$ by induction on $`n`$, starting from the single wait time and convolving in one more exponential at each step.

**Theorem.** For independent $`X_1, \dots, X_n \sim \text{Exponential}(\lambda)`$, the arrival time $`T_n = X_1 + \cdots + X_n`$ has density
```math
f_{T_n}(t) = \frac{\lambda^n t^{n-1}}{(n-1)!}\, e^{-\lambda t}, \qquad t \ge 0,
```
and $`f_{T_n}(t) = 0`$ for $`t < 0`$. This is the *Erlang density* with shape $`n`$ and rate $`\lambda`$.

**Note:** The Erlang density is the special case of the *Gamma distribution* in which the shape parameter is a positive integer. The Gamma distribution allows any real shape $`\alpha > 0`$, with $`(n-1)!`$ replaced by the Gamma function $`\Gamma(\alpha)`$ (which satisfies $`\Gamma(n) = (n-1)!`$); the Erlang is the name for the integer-shape members of that family, which are exactly the distributions of sums of independent exponentials.

**Proof.** The base case $`n = 1`$ is the exponential density itself: $`f_{T_1}(t) = f_X(t) = \lambda e^{-\lambda t}`$, which matches the formula since $`t^0 = 1`$ and $`0! = 1`$.

Assume the formula holds for some $`n \ge 1`$. Applying the convolution recursion with $`f_X(t - s) = \lambda e^{-\lambda (t - s)}`$,
```math
f_{T_{n+1}}(t) = \int_0^t \frac{\lambda^n s^{n-1}}{(n-1)!}\, e^{-\lambda s} \cdot \lambda e^{-\lambda(t-s)}\, ds.
```
The two exponentials combine to $`e^{-\lambda s} e^{-\lambda(t - s)} = e^{-\lambda t}`$, which does not depend on $`s`$ and comes out of the integral along with the constants:
```math
f_{T_{n+1}}(t) = \frac{\lambda^{n+1}}{(n-1)!}\, e^{-\lambda t} \int_0^t s^{n-1}\, ds.
```
The remaining integral is $`\int_0^t s^{n-1}\, ds = t^n / n`$, so
```math
f_{T_{n+1}}(t) = \frac{\lambda^{n+1}}{(n-1)!}\, e^{-\lambda t} \cdot \frac{t^n}{n} = \frac{\lambda^{n+1} t^n}{n!}\, e^{-\lambda t}.
```
This is the claimed formula with $`n`$ replaced by $`n+1`$, completing the induction. $`\quad\blacksquare`$

**Interpretation:** The two factors of the Erlang density carry the same reading as in the general wait-time picture. The exponential $`e^{-\lambda t}`$ is the survival factor, the chance of avoiding events long enough to still be accumulating them at time $`t`$. The polynomial $`\lambda^n t^{n-1}/(n-1)!`$ grows with $`t`$, reflecting that more elapsed time allows more ways for $`n`$ events to have fallen into $`[0, t]`$. Their product peaks and then decays: the $`n`$-th event is unlikely very early, when too little time has passed for $`n`$ events to accumulate, and unlikely very late, when survival that long is improbable.

## The Count $`N(t)`$ is Poisson Distributed

The counting identity relates the count to the arrival times, and the previous section shows the arrival times have Erlang density. With these two facts we can show that the distribution of $`N(t)`$ is the Poisson distribution from Axiom 3 of the definition of Poisson process.

### The Erlang survival function

The counting identity says $`P(N(t) \ge n) = P(T_n \le t)`$, and the right side is the Erlang cumulative distribution. Rather than integrate the Erlang density directly, we establish its cumulative distribution in a form that exposes the Poisson terms immediately.

**Lemma.** For $`n \ge 1`$ and $`t \ge 0`$, the arrival time $`T_n`$ satisfies
```math
P(T_n > t) = \sum_{k=0}^{n-1} e^{-\lambda t}\, \frac{(\lambda t)^k}{k!}.
```

**Proof.** Fix $`t`$ and induct on $`n`$. For $`n = 1`$, the right side is the single term $`e^{-\lambda t}`$, and the left side is $`P(T_1 > t) = P(X_1 > t) = e^{-\lambda t}`$, the exponential survival function; the two agree.

Assume the identity for some $`n \ge 1`$. The arrival times are nested, $`T_{n+1} = T_n + X_{n+1} \ge T_n`$, and passing from $`\{T_{n+1} > t\}`$ to $`\{T_n > t\}`$ differs by the event that the $`n`$-th arrival has occurred but the $`(n+1)`$-th has not:
```math
P(T_{n+1} > t) = P(T_n > t) + P(T_n \le t < T_{n+1}).
```
The second term is the probability of exactly $`n`$ arrivals by time $`t`$, which we compute from the Erlang density of $`T_n`$ and one more wait. Conditioned on $`T_n = s \le t`$, the $`(n+1)`$-th arrival exceeds $`t`$ exactly when the next wait exceeds $`t - s`$, an event of probability $`e^{-\lambda(t-s)}`$. Integrating against the Erlang density,
```math
\begin{align*}
P(T_n \le t < T_{n+1}) &= \int_0^t \frac{\lambda^n s^{n-1}}{(n-1)!}\, e^{-\lambda s}\, e^{-\lambda(t-s)}\, ds \\
&= \frac{\lambda^n e^{-\lambda t}}{(n-1)!} \int_0^t s^{n-1}\, ds \\
&= e^{-\lambda t}\, \frac{(\lambda t)^n}{n!}.
\end{align*}
```
Adding this to the inductive hypothesis for $`P(T_n > t)`$,
```math
\begin{align*}
P(T_{n+1} > t) &= \sum_{k=0}^{n-1} e^{-\lambda t}\, \frac{(\lambda t)^k}{k!} + e^{-\lambda t}\, \frac{(\lambda t)^n}{n!} \\
&= \sum_{k=0}^{n} e^{-\lambda t}\, \frac{(\lambda t)^k}{k!},
\end{align*}
```
which is the identity with $`n`$ replaced by $`n+1`$. $`\quad\blacksquare`$

### The Poisson distribution of $`N(t)`$

The lemma expresses the arrival-time distribution as a sum of Poisson terms, and the count now follows by differencing.

**Theorem.** For the counting process built from independent $`\text{Exponential}(\lambda)`$ wait times, and for each fixed $`t \ge 0`$,
```math
P(N(t) = n) = e^{-\lambda t}\, \frac{(\lambda t)^n}{n!}, \qquad n = 0, 1, 2, \dots,
```
that is, $`N(t) \sim \text{Poisson}(\lambda t)`$.

**Proof.** The count equals $`n`$ exactly when at least $`n`$ arrivals have occurred but not at least $`n+1`$:
```math
P(N(t) = n) = P(N(t) \ge n) - P(N(t) \ge n + 1).
```
By the counting identity, $`P(N(t) \ge n) = P(T_n \le t) = 1 - P(T_n > t)`$ for $`n \ge 1`$, while for $`n = 0`$ we have $`P(N(t) \ge 0) = 1`$ since the count is never negative. In both cases the lemma applies to the survival terms. For $`n \ge 1`$,
```math
\begin{align*}
P(N(t) = n) &= \big(1 - P(T_n > t)\big) - \big(1 - P(T_{n+1} > t)\big) \\
&= P(T_{n+1} > t) - P(T_n > t).
\end{align*}
```
The two sums from the lemma differ by exactly their last term,
```math
\begin{align*}
P(T_{n+1} > t) - P(T_n > t) &= \sum_{k=0}^{n} e^{-\lambda t}\frac{(\lambda t)^k}{k!} - \sum_{k=0}^{n-1} e^{-\lambda t}\frac{(\lambda t)^k}{k!} \\
&= e^{-\lambda t}\frac{(\lambda t)^n}{n!}.
\end{align*}
```
For $`n = 0`$, $`P(N(t) = 0) = 1 - P(T_1 \le t) = P(T_1 > t) = e^{-\lambda t}`$, which is the same formula at $`n = 0`$. This establishes the claim for all $`n \ge 0`$. $`\quad\blacksquare`$


Thus $`N(t) \sim \text{Poisson}(\lambda t)`$, establishing the $`s = 0`$ case of the third axiom.

## Independent and Stationary Increments

The previous section found the distribution of $`N(t)`$, the count measured from the origin. The definition asks for two more things: that the increment $`N(t) - N(s)`$ over an arbitrary interval $`(s, t]`$ is $`\text{Poisson}(\lambda(t-s))`$, and that increments over disjoint intervals are independent. Both follow from a single observation — the process, viewed from a fixed time $`s`$, is again the same construction — so we develop that observation first.

### The process restarts at $`s`$

Fix $`s \ge 0`$ and define the *shifted process*
```math
N_s(u) := N(s + u) - N(s), \qquad u \ge 0,
```
which counts the arrivals in $`(s, s + u]`$, renumbered so that time and count both restart at zero. We claim $`N_s`$ is a Poisson process with rate $`\lambda`$, built from its own sequence of independent $`\text{Exponential}(\lambda)`$ gaps, and that these gaps are independent of everything that happened before $`s`$.

The gaps of the shifted process are the waits between arrivals after $`s`$. All but the first are original gaps $`X_{N(s)+2}, X_{N(s)+3}, \dots`$, the waits that fall entirely after $`s`$; these are members of the original independent sequence $`\{X_n\}`$, so they are independent $`\text{Exponential}(\lambda)`$ variables and are independent of the history up to $`s`$, which involves only earlier gaps.

The only gap needing attention is the first, the residual wait $`R := T_{N(s)+1} - s`$ from $`s`$ to the next arrival. The arrival $`T_{N(s)+1}`$ is pending at $`s`$, so its gap $`X_{N(s)+1}`$ satisfies $`X_{N(s)+1} > s - T_{N(s)}`$: an amount $`a := s - T_{N(s)}`$ has already elapsed. The residual is $`R = X_{N(s)+1} - a`$, and by memorylessness
```math
\begin{align*}
P(R > u \mid X_{N(s)+1} > a) &= P(X_{N(s)+1} - a > u \mid X_{N(s)+1} > a) \\
&= P(X_{N(s)+1} > u) \\
&= e^{-\lambda u}.
\end{align*}
```
The residual is $`\text{Exponential}(\lambda)`$ regardless of the elapsed amount $`a`$, hence independent of it and of the earlier history. So $`R`$ too is a fresh $`\text{Exponential}(\lambda)`$ gap, independent of the past.

The shifted process is therefore assembled from independent $`\text{Exponential}(\lambda)`$ gaps in exactly the way the original process was, and all of these gaps are independent of the history up to $`s`$. 

**Note:** The residual step is the one place the construction uses the exponential specifically, and it is why the Poisson process is built on the exponential rather than some other wait time. Memorylessness makes the partially-elapsed straddling gap indistinguishable from a fresh one; §4 showed the exponential is the only continuous wait time for which this holds. With any other gap distribution the elapsed time $`a`$ would survive in the residual, and the process would not restart cleanly.

**Note:** The claim that $`R`$ is *independent of the history* because its conditional distribution does not depend on that history is stated here at the level of elementary conditioning, which suffices for restarting at a fixed time $`s`$. Making it fully rigorous requires describing the history as a $`\sigma`$-algebra $`\mathcal{F}_s`$ (the *filtration* generated by the process up to $`s`$) and showing the conditional distribution of the post-$`s`$ process given $`\mathcal{F}_s`$ is that of a fresh Poisson process; the care is warranted because the index $`N(s)+1`$ of the straddling gap is itself random. Restarting at a *random* time $`S`$ — for instance the time of the $`n`$-th arrival — is the *strong Markov property* and requires $`S`$ to be a *stopping time* together with a genuine measure-theoretic argument. A reader wanting the complete development should look under filtrations, the Markov property, and the strong Markov property in a text on stochastic processes.

### Stationary increments

Because $`N(t) - N(s) = N_s(t - s)`$ and $`N_s`$ is constructed from independent $\text{Exponential}(\lambda)$ wait times, the previous section's result applies to $`N_s`$, and the increment has the distribution computed there:
```math
N(t) - N(s) = N_s(t - s) \sim \text{Poisson}\big(\lambda (t - s)\big).
```
The distribution depends only on the interval length $`t - s`$, and the origin count $`N(t) \sim \text{Poisson}(\lambda t)`$ is the case $`s = 0`$. This is the definition's distributional clause in full.

### Independent increments

Take division points $`0 \le t_0 < t_1 < \cdots < t_k`$ with increments $`N(t_1) - N(t_0), \dots, N(t_k) - N(t_{k-1})`$. The shifted process $`N_{t_1}`$ is independent of the history up to $`t_1`$; the first increment $`N(t_1) - N(t_0)`$ is determined by that history, while every later increment is determined by $`N_{t_1}`$. So the first increment is independent of all the others together.

The later increments are increments of $`N_{t_1}`$, itself constructed from independent $`\text{Exponential}(\lambda)`$ wait times, so the same argument applied to $`N_{t_1}`$ at its first division point separates the second increment from those after it. Repeating down the list gives joint independence of the whole collection. $`\quad\blacksquare`$

## Uniqueness, Rate, and Varying the Rate

The construction from exponential wait times satisfies every clause of the definition, which settles existence: a Poisson process exists. Three tasks remain. We confirm the clauses determine the process uniquely, so that "the" Poisson process is a well-defined object; we identify what the rate $`\lambda`$ measures, connecting it to the hazard of the wait times; and we lift the constant rate to a rate that varies in time, recovering the general definition from the homogeneous case.

### Uniqueness

Existence shows the definition is not vacuous. Uniqueness shows it is not ambiguous: that any two processes meeting the definition have the same distribution, so the axioms describe one process rather than several. What must be checked is that the axioms fix every finite-dimensional distribution — the joint distribution of $`\big(N(t_1), \dots, N(t_k)\big)`$ for any finite set of times — since a stochastic process is determined by these joint distributions.

Fix times $`0 = t_0 < t_1 < \cdots < t_k`$. The counts and the increments carry the same information: knowing the increments $`D_i := N(t_i) - N(t_{i-1})`$ is equivalent to knowing the counts, through $`N(t_i) = D_1 + \cdots + D_i`$ (using $`N(t_0) = N(0) = 0`$). It therefore suffices to determine the joint distribution of $`(D_1, \dots, D_k)`$. The definition does this directly: the increments are independent (clause 2), so their joint distribution is the product of their marginals, and each marginal is fixed by the distributional clause (clause 3),
```math
D_i = N(t_i) - N(t_{i-1}) \sim \text{Poisson}\big(\lambda(t_i - t_{i-1})\big).
```
The joint distribution is thus
```math
P(D_1 = n_1, \dots, D_k = n_k) = \prod_{i=1}^k e^{-\lambda(t_i - t_{i-1})}\frac{\big(\lambda(t_i - t_{i-1})\big)^{n_i}}{n_i!},
```
completely determined by $`\lambda`$ and the chosen times. Any two processes satisfying the definition produce this same joint distribution at every finite collection of times, so they are equal in distribution. The Poisson process of rate $`\lambda`$ is unique.

### What the rate measures

The parameter $`\lambda`$ entered twice, and the two entries are worth holding side by side. In the wait-time half of the chapter, $`\lambda`$ was the *hazard* of each exponential gap: the instantaneous rate of the next event given that it has not yet occurred. In the process, $`\lambda`$ turns out also to be the mean count per unit time, since $`E[N(t)] = \lambda t`$ gives
```math
\frac{E[N(t)]}{t} = \lambda.
```
These are two readings of one number, and their agreement is a feature of the homogeneous case rather than a coincidence. The hazard is a conditional rate — events per unit time among gaps still in progress — while the mean count per unit time is an unconditional rate, averaged over the whole process. They coincide because memorylessness erases the distinction: a gap in progress is distributed like a fresh one, so conditioning on survival changes nothing, and the conditional and unconditional rates are equal. The single word *rate* serves both because, for the exponential, both are $`\lambda`$.

This is also why the distinction cannot be ignored once the rate varies. When the hazard changes with time, a gap in progress is no longer like a fresh one, and the instantaneous rate at time $`t`$ need not equal the count accumulated per unit time. Keeping the two readings straight is what makes the inhomogeneous process well-posed.

### Varying the rate

The construction produces only the homogeneous process: stacking i.i.d. exponential gaps gives a single rate $`\lambda`$, the same at every time. A rate that varies with time breaks the "identically distributed" assumption and with it the clean sum-of-exponentials structure, so we do not reach the inhomogeneous process by rebuilding the stack. We reach it by bending time.

**Def:** Let $`\lambda : [0, \infty) \to [0, \infty)`$ be a rate function. Its accumulated total, the *cumulative rate function*, is defined by
```math
\Lambda(t) := \int_0^t \lambda(u)\, du.
```

$`\Lambda`$ is nondecreasing, and it measures elapsed time not in ordinary units but in accumulated rate — the same reshaping of the time axis by an integrated rate that appeared as the cumulative hazard in the wait-time chapter, where $H_X(x) = \int_0^x h_X(w)\,dw$ turned any lifetime into a standard one. Here $`\Lambda`$ plays the same role for a whole process.

**Def:** Let $`M`$ be a homogeneous Poisson process of rate $`1`$. The *inhomogeneous Poisson process* with rate function $`\lambda`$ is
```math
N(t) := M\big(\Lambda(t)\big),
```
the rate-$`1`$ process read on the warped clock $`\Lambda`$.

Running the rate-$`1`$ process on the clock $`\Lambda`$ speeds it up where $`\lambda`$ is large and slows it where $`\lambda`$ is small, so that the local rate of $`N`$ at time $`t`$ is $`\lambda(t)`$. The axioms transfer from $`M`$ with almost no work, because $`\Lambda`$ is deterministic and monotone.

**Theorem.** The process $`N(t) = M(\Lambda(t))`$ satisfies the definition of a Poisson process with rate function $`\lambda(t)`$: it has $`N(0) = 0`$, independent increments, and $`N(t) - N(s) \sim \text{Poisson}\big(\int_s^t \lambda(u)\,du\big)`$.

**Proof.** Since $`\Lambda(0) = 0`$ and $`M(0) = 0`$, we have $`N(0) = 0`$.

For increments, an interval $`(s, t]`$ in ordinary time maps to the interval $`(\Lambda(s), \Lambda(t)]`$ in clock time, and
```math
N(t) - N(s) = M(\Lambda(t)) - M(\Lambda(s)).
```
Because $`\Lambda`$ is nondecreasing, disjoint intervals $`(s_i, t_i]`$ map to disjoint clock intervals $`(\Lambda(s_i), \Lambda(t_i)]`$. The increments of $`N`$ over the original intervals are exactly the increments of $`M`$ over the images, and $`M`$ has independent increments over disjoint intervals, so the increments of $`N`$ are independent.

For the distribution, $`M`$ has rate $`1`$, so by the homogeneous result its increment over a clock interval of length $`\Lambda(t) - \Lambda(s)`$ is
```math
M(\Lambda(t)) - M(\Lambda(s)) \sim \text{Poisson}\big(\Lambda(t) - \Lambda(s)\big) = \text{Poisson}\!\left(\int_s^t \lambda(u)\, du\right),
```
using $`\Lambda(t) - \Lambda(s) = \int_s^t \lambda(u)\,du`$. This is the required increment distribution. $`\quad\blacksquare`$

The homogeneous process is the special case $`\lambda(t) \equiv \lambda`$, where $`\Lambda(t) = \lambda t`$ and the warped clock is a linear rescaling: $`N(t) = M(\lambda t)`$ runs the unit-rate process at constant speed $`\lambda`$. The increment distribution collapses to $`\text{Poisson}(\lambda(t - s))`$, and the general definition returns to the one we proved from stacked exponentials.
## The Infinitesimal View

The chapter has told its whole story through gaps: wait times stacked into arrivals, arrivals counted into $`N(t)`$. A complementary view looks at the process through a vanishingly small window of time rather than across a full waiting interval. It recovers two facts the gap picture leaves implicit — that events never occur simultaneously, and that the Poisson count is what one expects from many independent rare chances — and it exposes a second, equivalent way to define the process.

Everything here reads off the increment distribution already in hand. We work in the homogeneous case for clarity; the inhomogeneous case follows by replacing $`\lambda h`$ with $`\int_t^{t+h}\lambda(u)\,du \approx \lambda(t) h`$ throughout. Consider the count over a short interval $`(t, t+h]`$, which is $`\text{Poisson}(\lambda h)`$. Expanding the Poisson probabilities for small $`h`$,
```math
\begin{aligned}
P\big(N(t+h) - N(t) = 0\big) &= e^{-\lambda h} = 1 - \lambda h + o(h), & &\text{(1)}\\
P\big(N(t+h) - N(t) = 1\big) &= \lambda h\, e^{-\lambda h} = \lambda h + o(h), & &\text{(2)}\\
P\big(N(t+h) - N(t) \ge 2\big) &= o(h). & &\text{(3)}
\end{aligned}
```
In a short window the process almost always does nothing, occasionally records a single event with probability proportional to the window's width, and essentially never records more than one. Statements (1)–(3) are the infinitesimal view.

### Events do not coincide

Statement (3) says the probability of two or more events in a window of width $`h`$ is negligible compared to $`h`$ itself. Shrinking the window, the chance of a coincidence vanishes faster than the window does, so in the limit no instant carries more than one event. The process increases one step at a time.

This agrees with what the construction already guaranteed. The arrival times are strictly increasing when the wait times are continuous, since each gap is positive with probability one, so two arrivals cannot share an instant. The gap picture forbids coincidences because positive gaps separate the arrivals; the infinitesimal picture forbids them because $`P\big(N(t+h) - N(t) \ge 2\big) = o(h)`$. The two views reach the same conclusion from opposite ends.

### Poisson from many rare chances

Statements (1) and (2) give a second route to the Poisson count, one that explains the distribution's name and matches the classical intuition. Partition an interval $`(s, t]`$ into $`m`$ subintervals of equal width $`h = (t - s)/m`$. Over each subinterval the process, by the expansion above, records one event with probability $`p = \lambda h + o(h)`$ and no event otherwise, the two-or-more case being negligible. Each subinterval is then nearly a Bernoulli trial: event or no event. By independent increments, the trials over disjoint subintervals are independent.

The count over $`(s, t]`$ is the total number of successes across the $`m`$ trials, so it is approximately $`\text{Binomial}(m, p)`$, with expected value
```math
m p = m\big(\lambda h + o(h)\big) = \lambda (t - s) + o(1).
```
As $`m \to \infty`$ the subintervals shrink, the per-trial probability $`p \to 0`$, and the number of trials $`m \to \infty`$, while the expected count $`m p`$ stays fixed at $`\lambda(t - s)`$. This is the regime of the following classical result.

**Theorem (Poisson limit theorem).** If $`Y_m \sim \text{Binomial}(m, p_m)`$ with $`m p_m \to \mu`$ as $`m \to \infty`$, then $`Y_m \to \text{Poisson}(\mu)`$ in distribution:
```math
\lim_{m\to\infty} P(Y_m = k) = e^{-\mu}\frac{\mu^k}{k!}, \qquad k = 0, 1, 2, \dots
```

Applying it with $`\mu = \lambda(t-s)`$ recovers $`N(t) - N(s) \sim \text{Poisson}(\lambda(t-s))`$, the increment distribution we derived earlier from stacked exponentials. The two derivations show different things. The stacked-exponential proof shows the Poisson count is forced by memoryless gaps; the rare-chance limit shows it is the count of many independent, individually unlikely opportunities for an event — the sense in which the Poisson distribution is a law of rare events.

### The infinitesimal characterization

The expansion that opened this section is itself a second definition. A process satisfying $`N(0) = 0`$, independent increments, and the infinitesimal conditions
```math
\begin{aligned}
P\big(N(t+h) - N(t) = 1\big) &= \lambda(t)\, h + o(h), \\
P\big(N(t+h) - N(t) \ge 2\big) &= o(h),
\end{aligned}
```
(with the zero-event probability determined as the remainder) is a Poisson process in the sense of our original definition, and conversely. The two definitions are equivalent.

**Note:** We have shown one direction: our process, defined through Poisson-distributed increments, satisfies the infinitesimal conditions, since they are the small-$`h`$ expansion of those increments. The converse — that a process meeting the infinitesimal conditions has Poisson-distributed increments — is proved by converting the infinitesimal conditions into a differential equation for the increment probabilities $`P(N(t+h) - N(s) = k)`$ and solving it. The result is a system whose solution is the Poisson distribution, much as the constant-hazard differential equation earlier in the chapter had the exponential as its solution. A reader wanting the full argument should look for the derivation of the Poisson process from its infinitesimal rates in a text on stochastic processes.

The infinitesimal conditions are often taken as the *starting* definition of the Poisson process, with the Poisson distribution of the counts derived from them. Our development ran the other way, building the process from memoryless wait times and proving the counts are Poisson.
