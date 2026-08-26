# 1.0 Parameterizing and Modeling Transmission

## Simple Poisson Process

Assume a person has intrinsic infectiousness of 1 person per day for time $`[0, 2]`$.

By this, we mean that the conditional hazard is 1 transmission attempt per day:

```math
\displaystyle \lim_{\varepsilon \to 0} \frac{P(\text{transmit in } [t, t + \varepsilon] \mid \text{have not transmitted by time } t)}{\varepsilon}=1.
```

For one fixed susceptible contact, the first successful transmission time is a waiting-time random variable. For infection attempts across all possible contacts, we instead use a Poisson point process, which can have multiple events.


### Q: How would you simulate the probability and time of infection for a single susceptible contact?

**Answer:** Sample $`\Delta t \sim \text{Exp}(1)`$. If $`\Delta t \lt  2`$, then the time between the infector becoming infected and this first infection is $`\Delta t`$. The person is infected at time $`\Delta t`$ after the infector was infected. Otherwise, the person is not infected.

### Q: How would you simulate a Poisson process with this same rate function?

```math
r(t) = \begin{cases}
1, & t \in [0, 2] \\
0, & \text{else}
\end{cases}
```

**Option 1:** First compute the number of events, then distribute in time.
Sample from a $`\text{Poisson}(2) \to n`$. Then sample $`n`$ times uniformly in $`[0, 2]`$ to get the times for the events.

*Intuition:* Case $`n = 1`$: We are assuming $`\Delta t_1 \lt  2`$ AND $`\Delta t_1 + \Delta t_2 > 2`$.

**Option 2:** Iteratively sample inter event times $`\Delta t_i \sim \text{Exp}(1)`$. Then the times are $`\Delta t_1, \Delta t_1 + \Delta t_2, \Delta t_1 + \Delta t_2 + \Delta t_3,\dots`$ up to when they exceed 2 (past which, our rate function says that no infection attempts are expected with a rate of 0).

**Option 3**: Time scaling. We can sample the inter event times using the cumulative rate function $\Lambda(t)$. For rate function $r(t)$ (equivalent to the hazard function in this case)

```math
r(t) = \begin{cases}
1, & t \in [0, 2] \\
0, & \text{else}
\end{cases}
```

the cumulative rate function $\Lambda(t)$ is

```math
\Lambda(t) = \begin{cases}
t, & t \in [0, 2] \\
2, & else
\end{cases}
```

A true inverse function of $\Lambda(t)$ does not exist on $[0,\infty)$ because $\Lambda(t)=2$ for every $t\geq2$. However, restricted to cumulative values $y\in[0,2]$, we can use $d(y)=y$. On this relevant range, $\Lambda(d(y))=y$; and for time values $t\in[0,2]$, $d(\Lambda(t))=t$.

Recall that inter event distances in the cumulative space can be sampled as $\Delta y_i \sim \text{Exp(1)}$. Iteratively sample $\Delta y_i$. Then the points in cumulative space are $\Delta y_1, \Delta y_1 + \Delta y_2, \Delta y_1 + \Delta y_2 + \Delta y_3, ...$, up to when they exceed 2 (the maximum value of $\Lambda(t)$, past which there is no valid inverse function). We'll call these $y_i$.

Use $d(y)$ to map $y_i$ from the cumulative space back to the time space, i.e. $t_i = d(y_i)$. These are the forecasted times for infection events relative to the infection time of the infector. The inter event times can now be calculated as $\Delta t_i = t_i - t_{i-1}$.

### Observation:
If we simulate from a Poisson process of rate $`r(t)`$, then the probability of having one event and the time of that event matches the intrinsic infectious process.


## Rejection Sampling and Time Scaling

Consider how to handle changes:
Suppose at a random time $`t \in [0, 2]`$, this person will wear a perfect facemask.

### Q: How does this change intrinsic transmissibility, and how do we simulate?

**Option 1:** Rejection sample. You sample according to Option 2. Then at the time of the purported event, if it is after the facemask time $`t`$, you reject it (ignore it).

What if the facemask isn’t perfect? Say it only prevents transmission 30% of the time.

**Answer:** You still schedule it normally, and then reject with probability 0.3.

*Observation:*
This 30% is not the reduction in the probability of infecting a single contact. It is the reduction in the probability of a single infection attempt succeeding.

This is called **Poisson thinning**: retaining each proposed event with probability $`0.7`$ changes the attempt rate from $`r(t)`$ to $`0.7r(t)`$.


**Option 2**: Time scaling. Rather than sampling inter event times iteratively directly as $x_i\sim \text{Exp(1)}$ and summing $x_i$ to get event times, let's sample the inter event distances in the cumulative space as $\Delta y_i \sim \text{Exp(1)}$.

Then calculate $y_i = \Sigma_{i' = 1}^{i} \Delta y_{i'}$

Each $y_i$ can be inverted with $d(t)$ to get $t_i = d(y_i)$.

Then the inter event times can be calculated as $\Delta t_i = t_i - t_{i-1}$, where $t_0 = 0$ for the time of the zeroeth event.

When an event causes the rate function to change after some time $`t^{*}`$, rather than rejecting events after $`t^{*}`$ with probability $`30\%`$, instead we re-evaluate what each value $y_i$ in the cumulative space maps to in the time space for events $t_i$ scheduled to happen after $`t^{*}`$. Since the rate function has changed after $`t^{*}`$, for those events we calculate $t_i$ with the inverse of the new cumulative rate function, $`t_i = d^{*}(y_i)`$, and calculate the inter event times as $\Delta t_i = t_{i} - t_{i-1}$. Reschedule events with the new event times.

*Observation:*
If the person took an antiviral that changes their new infectiousness, you could assume this is the person’s new intrinsic infectiousness. Sometimes, you may *have* to:
Imagine viral rebound. We have to reschedule if it *increases* future transmissibility. In this case, time scaling and forecasting events through the cumulative rate function is a straight forward approach.

Whenever rescheduling is needed, start from the time of the change, not from the original scheduling time. Suppose the rate changes at time $`t^{*}`$ before the next proposed attempt. Because that proposal is still in the future, we know that no **proposed forecast event** occurred before $`t^{*}`$; it is not enough merely to know that no actual transmission occurred. In cumulative-rate space, the remaining distance to the next proposed attempt is still $`\text{Exp}(1)`$ by memorylessness. Therefore, cancel the pending proposal, reset the cumulative-rate clock at $`t^{*}`$, use the new upper-bound forecast after $`t^{*}`$, and sample a new next attempt from there. The same logic applies after a proposed attempt is rejected by thinning: future Poisson-process increments are independent of the rejected attempt, so generate the next proposal starting from the rejection time. A rate increase requires this rescheduling because the old forecast may no longer be an upper bound on the actual rate.

## Distributing Infection Hazard

Using the example above of a person who is uniformly infectious on $`[0, 2]`$,
suppose an infected person lives in a house with 2 other people, and they stay in that house for all time during the course of the infection.

### How should we model the infection process of those two people in terms of an inhomogeneous Poisson process?

**Possibility 1:** The hazard of infecting each person is equal and half of the intrinsic infectiousness. The rate of the total infection process is, in this case, just the intrinsic infectiousness (halved and then added for each person).

**Possibility 2 :** The hazard of infecting each person is equal to the intrinsic infectiousness. The rate process of the total infection process is twice the intrinsic infectiousness.

In the real world we observe something between these two extremes.

### Example from History:
We assumed there was a magic parameter $`\alpha`$ and that in a setting with $`n`$ people (total, including the source of infection), the rate of the total infectiousness process was:

```math
\text{Intrinsic Infectiousness} \times (n - 1)^\alpha \quad \text{where} \quad 0 \leq \alpha \leq 1,
```

which interpolates between the two options.


## Modeling Transmission with Time-Varying Factors

Suppose this person, who is in a house with 2 people, is going to wear a facemask at random times throughout the day.

### Claim:
Let $`r(t)`$ be the intrinsic infectiousness for this person in the absence of any interventions. You schedule infections according to $`2^\alpha \cdot r(t)`$ and reject them with probability 0.3 if the person is wearing a facemask at that time and otherwise accept. This works because $`r(t) \mapsto 2^\alpha \cdot r(t)`$ is linear in $`r(t)`$.

However, total infectiousness in this setting could be some function $`f_2(r(t))`$ for some arbitrary monotonic positive function. In this case, you would accept with probability:

```math
\frac{f_2(0.7 \cdot r(t))}{f_2(r(t))}
```

### Linearity is Convenient:
Total infectiousness can be modeled as:

```math
\text{Total infectiousness} = \left( \text{fraction of time in home} \cdot \text{home multiplier} + \text{fraction of time in work} \cdot \text{work multiplier} \right) \cdot r(t)
```

Suppose a person spends half of their time at home and half at work. Home has 3 people, work has 9 people. Assume $`\alpha = 1`$ as above.
The total infectiousness would be
$`\left( 0.5 \cdot 2 + 0.5 \cdot 8 \right) \cdot r(t)`$.

### How Would We Handle Telework?
From the perspective of the source, you adjust the fraction of time.

**Proposal:** Schedule infection attempts at a rate of
$`(0.5 \cdot 2 + 0.5 \cdot 8) \cdot r(t)`$
and record the expected total infectiousness at that time in the plan. Call that $`R`$.
At the time of the attempt, compute the person’s total infectiousness in each setting. Let $`S`$ be the sum of these infectiousnesses.
You accept the proposed attempt with probability
$`S/R`$.
Then, you choose the setting proportional to their current total infectiousness. After that, you pick the contact in that setting.

Assume a person is either at home or at work at any point in time, but you don’t know the schedule.
At home, their total infectiousness would be $`2 \cdot r(t)`$. At work, it would be $`8 \cdot r(t)`$.

**Proposal:** Schedule based on $`8 \cdot r(t)`$. Proceed as above.

**Why use the maximum?** This is a thinning algorithm, so the scheduled forecast rate $`R(t)`$ must be an upper bound on the actual total infectiousness $`S(t)`$ at every possible attempt time. Only then is the acceptance probability $`S(t)/R(t)`$ between $`0`$ and $`1`$. Here the person could be at work, where $`S(t)=8r(t)`$, so $`R(t)=8r(t)`$ is a valid forecast. The expected rate $`(0.5\cdot2+0.5\cdot8)r(t)=5r(t)`$ is not valid when the person's location is unknown: if they are at work, it would require accepting with probability $`8r(t)/5r(t)=8/5>1`$.

More importantly, a forecast below the actual rate can omit attempts that should have occurred. Rejecting excess proposed attempts is safe; trying to create missing attempts later is not. If new information can increase the future actual rate above the current forecast, cancel the pending attempt and reschedule using a new upper-bound forecast.

## Algorithm
## Sketch of Transmission Model Implementation

### Concepts

1. **Transmission settings**
   Conceptual locations where transmission may occur, e.g.
   1. Home, School, Work, Community

2. **Itineraries**
   A potentially time-varying proportion of time allocation for people in various transmission settings, e.g.
   1. 40% Home; 40% Work; 20% Community all time
   2. 100% Work 9a-5pm; 80% Home, 20% Community otherwise

3. **Infectiousness functions**
   The rate function of a Poisson process that represents infections or infection attempts, potentially in a setting or across multiple settings.

### Steps

1. **Person is infected**
2. **Data is recorded** about the person that will govern their future infectiousness, symptoms, and behavior, potentially scheduling plans.
3. **Infection propagation loop:**
   1. Determine a forecast of the future infectiousness of a person based on currently-known information. This should be some infectiousness function that is at least as large as the true total infectiousness of the person at all future points in time.
   2. Use this forecast to determine the next transmission attempt time and the forecasted infectiousness at that time. If this time is infinite, the person can be assumed to no longer be infectious, and the infection propagation loop ends.
   3. Schedule a transmission attempt for this time and record the forecasted infectiousness at that attempt time.
   4. At the time of the transmission attempt, compute the current infectiousness for the source person across each setting as well as its total.
   5. Accept the transmission attempt with the probability of the ratio of the current total infectiousness divided by the forecasted infectiousness. If accepted:
      1. Randomly select the setting of transmission proportionally to the current infectiousness in each setting.
      2. Randomly select a contact for transmission in that setting.
      3. Compute the current susceptibility of the contact, and infect them with this probability.
   6. Go to step 1.

4. If an event occurs that will increase the forecasted future rate of infection between the time of this event and the next transmission attempt time, then the currently-scheduled infection attempt must be cancelled (if it exists) and recomputed from step 1 above with the new forecast.

### Example to Work Through:

**Two transmission settings:** Home, Work

All people spend time in Home and Work according to some time schedule, potentially random. Suppose

```math
\text{intrinsic infectiousness} =
\begin{cases}
1, & 0 \leq t \lt  1 \\
2, & 1 \leq t \lt  2 \\
0, & t \geq 2
\end{cases}
```

For a specific person who was just infected, assume:
- Their Home setting has up to 3 people in it (including themselves).
- Their Work setting has up to 5 people in it (including themselves).

Assume infection risk applies independently to people in all settings ($`\alpha = 1`$ as above).

**Loop:**

The intrinsic infectiousness is at most 2, so the total infectiousness is at most 4 at Home and at most 8 at Work. We can use either of the following upper-bound forecasts.

**Case A: Constant upper-bound forecast.** Use $`R(t)=8`$ for $`0\leq t\leq2`$.

1. Sample $`X\sim\text{Exp}(1)`$ and set $`t=X/8`$.
2. If $`t>2`$, there is no proposed attempt; stop.
3. Otherwise, propose an attempt at time $`t`$ with forecast rate $`R(t)=8`$.

**Case B: Piecewise upper-bound forecast.** Use $`R(t)=4`$ for $`0\leq t<1`$ and $`R(t)=8`$ for $`1\leq t<2`$.

1. Sample a cumulative-rate value $`X\sim\text{Exp}(1)`$. We will choose a time $`t`$ such that $`\Lambda(t)=X`$; equivalently, $`t=\Lambda^{-1}(X)`$.
2. The cumulative forecast rate is $`\Lambda(t)=4t`$ for $`0\leq t\leq1`$ and $`\Lambda(t)=4+8(t-1)=8t-4`$ for $`1\leq t\leq2`$. Its total is $`\Lambda(2)=12`$.
3. If $`X>12`$, there is no proposed attempt; stop. Otherwise, map back to time:

```math
t=d(X)=
\begin{cases}
X/4, & 0\leq X\leq4,\\
(X+4)/8, & 4<X\leq12.
\end{cases}
```

4. Propose an attempt at time $`t`$ with forecast rate $`R(t)`$.

For every finite proposed time in Case B, $`X=\Lambda(t)`$: $`X`$ is the cumulative-rate value, and $`t=\Lambda^{-1}(X)`$ is the corresponding time. If $`X>12`$, no time within the infectious period corresponds to that cumulative-rate value.

**Shared steps after either proposal:**

1. At the proposed time, calculate the actual total infectiousness $`S(t)`$ and accept the attempt with probability $`S(t)/R(t)`$.
2. For example, suppose the proposed time is $`t=0.5`$. If the person is at Home with only one other person present, then $`S(t)=1`$. The acceptance probability is $`1/8`$ in Case A and $`1/4`$ in Case B. If they are at Work with three other people present, then $`S(t)=3`$, giving acceptance probabilities $`3/8`$ and $`3/4`$, respectively.
3. If the attempt is accepted, choose the setting and contact, then infect the contact if they are susceptible.

### Q: What if people can wear a facemask?

Treat a facemask as a time-varying reduction in the actual infectiousness. For example, if it reduces transmission by a fraction $`e`$ whenever it is worn, replace $`S(t)`$ by

```math
S_{\mathrm{mask}}(t)=(1-e)S(t)
```

in the shared acceptance step. The proposed rate $`R(t)`$ need not change when the mask can only reduce infectiousness; accept with probability $`S_{\mathrm{mask}}(t)/R(t)`$. This is again Poisson thinning. If new information could *increase* the actual rate above the current forecast $`R(t)`$, the pending proposal must instead be cancelled and resampled with a larger upper-bound forecast.

If mask effectiveness differs by setting or contact, first calculate the adjusted infectiousness in each setting. Accept with probability

```math
\frac{\text{sum of adjusted infectiousnesses across settings}}{R(t)},
```

then choose the setting proportionally to its adjusted infectiousness and choose a contact within that setting.

---

### Negative Binomial Branching Process Model

Suppose $`f(t)`$ is a generation-interval density, so $`\int_0^\infty f(t)\,dt=1`$. We want each infected person to have a negative-binomial number of secondary infections while the timing of each secondary infection follows $`f`$.

**Intrinsic infectiousness:** For a person with individual infectiousness multiplier $`A`$, use the rate function

```math
r(t)=A f(t).
```

Conditional on $`A`$, the total rate integrates to

```math
\int_0^\infty A f(t)\,dt=A.
```

Therefore, conditional on a given person's infectiousness multiplier $`A`$, their number of secondary infections $`K`$ is Poisson distributed:

```math
K\mid A\sim\text{Poisson}(A).
```

Thus $`E[K\mid A]=A`$: for a given person, $`A`$ is their expected number of secondary infections. Averaging over people gives $`E[K]=E[A]`$, the population-average expected number of secondary infections.

To obtain overdispersion, sample a different $`A`$ for each infected person from a Gamma distribution. If

```math
A\sim\text{Gamma}\left(k,\;\text{rate}=\frac{k}{R_0}\right),
```

then $`E[A]=R_0`$ and the marginal offspring distribution is negative binomial with mean $`R_0`$ and dispersion $`k`$:

```math
K\sim\text{NegBin}(\text{mean}=R_0,\;\text{dispersion}=k).
```

To simulate a person's transmissions, either simulate a Poisson point process with rate $`A f(t)`$, or equivalently:

1. Sample $`A`$ from the Gamma distribution.
2. Sample $`K\sim\text{Poisson}(A)`$.
3. Sample $`K`$ independent transmission times from density $`f(t)`$.

### Worked Example: Negative-Binomial Infectiousness with a Facemask

Suppose a person's individual multiplier is $`A`$ and their baseline rate is $`A f(t)`$. They wear a facemask from time $`0.5`$ through time $`1.2`$, and the mask retains a fraction $`m=0.7`$ of the baseline transmission rate. Define

```math
m(t)=
\begin{cases}
0.7, & 0.5\leq t\leq1.2,\\
1, & \text{otherwise}.
\end{cases}
```

The actual rate is

```math
r_{\mathrm{actual}}(t)=A f(t)m(t).
```

#### Direct sampling when it is easy

Let

```math
C=\int_0^\infty f(t)m(t)\,dt.
```

Conditional on $`A`$, the expected number of secondary infections is now $`AC`$, so

```math
K\mid A\sim\text{Poisson}(AC).
```

Conditional on $`K`$, the transmission times are independent draws from the adjusted density

```math
g(t)=\frac{f(t)m(t)}{C}.
```

If $`C`$ is easy to calculate and there is a direct sampler for $`g`$, this is the simplest method: sample $`A`$, sample $`K`$, then sample $`K`$ times from $`g`$. In the special case where the mask is worn for the whole infection, $`m(t)=0.7`$, so $`C=0.7`$ and $`g(t)=f(t)`$: only the expected count changes.

#### Forecasting and thinning when the rate can changes

For a more realistic simulation, the mask schedule may be unknown until later, may change repeatedly, or may depend on a person's location. Instead use the baseline rate as an upper-bound forecast:

```math
R(t)=A f(t).
```

Generate proposed attempts from this forecast. Its cumulative forecast rate is

```math
\Lambda_R(t)=\int_0^t A f(s)\,ds=A F(t).
```

To generate proposed times by time scaling, sample independent $`E_i\sim\text{Exp}(1)`$ and form $`Y_i=E_1+\cdots+E_i`$. If $`Y_i>A`$, stop: the total cumulative forecast rate is $`A`$. Otherwise, map the cumulative-rate value back to time:

```math
t_i=\Lambda_R^{-1}(Y_i)=F^{-1}\left(\frac{Y_i}{A}\right).
```

At each proposed time $`t`$, check whether the mask is actually being worn and accept the proposal with probability $`m(t)`$. The accepted attempts then have rate

```math
R(t)m(t)=A f(t)m(t)=r_{\mathrm{actual}}(t).
```

This is valid because $`0\leq m(t)\leq1`$: the forecast never falls below the actual rate. A later mask change that only reduces transmission needs no rescheduling; it merely changes the acceptance probability of future proposed attempts. A later change that can increase the rate above $`R(t)`$ requires a new, larger forecast and rescheduling.

If $`g`$ is easy to sample, direct sampling is shorter and more efficient because it generates only actual events. Forecasting/time scaling is useful when $`g`$ has no convenient sampler or must be repeatedly recomputed as behavior, interventions, contacts, or locations change. It only requires a usable cumulative forecast rate and a way to invert it (often numerically), while thinning safely discards proposed attempts that no longer apply.
