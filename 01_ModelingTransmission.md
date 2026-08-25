# 1.0 Parameterizing and Modeling Transmission

## Simple Poisson Process

Assume a person has intrinsic infectiousness of 2 people per day for time $`[0, 1]`$.

By definition, this means that if $`t \in [0,1)`$, and $`\Delta t`$ is very small, then the probability of transmitting in $`[t, t+\Delta t]`$ is approximately $`2 \Delta t`$.  More precisely the error is making this approximation is small compared to $`\Delta t`$ for small enough $`\Delta t`$:

$`\lim_{\Delta t \to 0} \frac{P(\text transmit in [t+\Delta t])}{\Delta t} = 2`$


### Q: How would you simulate the probability and time of infection for a single susceptible contact (for which we care only about the first event)?

**Answer:** Sample $`\Delta t \sim \text{Exp}(2)`$. If $`\Delta t \lt  1`$, then the waiting time between the infector becoming infected and this first infection is $`\Delta t`$. The person is infected at time $`\Delta t`$ after the infector was infected. Otherwise, the person is not infected.

### Q: How would you simulate a Poisson process (for which multiple events can occur) with this same rate function?

```math
r(t) = \begin{cases}
2, & t \in [0, 1] \\
0, & \text{else}
\end{cases}
```

**Option 1:** First compute the number of events, then distribute in time.  

Given a Poisson process of rate $2$ for $1$ units of time, the number that occur have a Poisson distribution. 

Sample the number of events $`n`$ from a $`\text{Poisson}(2)`$ distribution. Then sample $`n`$ times uniformly in $`[0, 1]`$ to get the event times.

*Intuition:* Case $`n = 1`$: We are assuming $`\Delta t_1 \lt  1`$ AND $`\Delta t_1 + \Delta t_2 > 1`$.

**Option 2:** Iteratively sample inter event times $`\Delta t_i \sim \text{Exp}(2)`$. Then the times are $t_1=`\Delta t_1`$,  $`t_2=t_1 + \Delta t_2`$,  $`t_3= t_2 + \Delta t_3$, $\ldots`$ up to when they exceed 1 (past which, our rate function says that no infection attempts are expected with a rate of 0).

**Option 3**: Time scaling. We can sample the inter event times using the cumulative rate function $c(t)$. For rate function $r(t)$ (equivalent to the hazard function in this case)

```math
r(t) = \begin{cases}
2, & t \in [0, 1] \\
0, & \text{else}
\end{cases}
```

the cumulative rate function $c(t)$ is

```math
c(t) = \begin{cases}
2t, & t \in [0, 1] \\
2, & else
\end{cases}
```

A true inverse function of $c(t)$ doesn't exist $\forall t \in [0, \infty)$, but over the interval $[0, 1]$ the inverse can be defined as $d(y)$ = y/2. (Why? If $d(y) = y/2$, then $d(c(t)) = c(t)/2 = 2t/2=t$ and similarly $c(d(t)) = 2d(t) = 2t/2=t$.)

Recall that inter event distances in the cumulative space can be sampled as $\Delta y_i \sim \text{Exp(1)}$. Iteratively sample $\Delta y_i$. Then the points in cumulative space are $\Delta y_1, \Delta y_1 + \Delta y_2, \Delta y_1 + \Delta y_2 + \Delta y_3, ...$, up to when they exceed 2 (the maximum value of $c(t)$, past which there is no valid inverse function). We'll call these $y_i$. 

Use $d(t)$ to map $y_i$ from the cumulative space back to the time space, i.e. $t_i = d(y_i)$. These are the forecasted times for infection events relative to the infection time of the infector. The inter event times can now be calculated as $\Delta t_i = t_i - t_{i-1}$.

### Observation:
If we simulate from a Poisson process of rate $`r(t)`$, then the probability of having one event and the time of that event matches the intrinsic infectious process.


## Rejection Sampling and Time Scaling

Consider how to handle changes:
Suppose at some time $`t^{\text{*}} \in [0, 1]`$, this person will wear a perfect facemask.

### Q: How does this change intrinsic transmissibility, and how do we simulate?

**Option 1:** Rejection sample. You sample as before (say, according to Option 2 above). Then at the time of the purported event, if it is after the facemask time $`t^{\text{*}}`$, you reject it (ignore it).

What if the facemask isn’t perfect? Say it only prevents transmission 30% of the time.

**Answer:** You still schedule it normally, and then reject with probability 0.3.

*Observation:*
This 30% is not the reduction in the probability of infecting a single contact. It is the reduction in the probability of a single infection attempt succeeding.


**Option 2**: Time scaling. Rather than sampling inter event times iteratively directly as $\Delta t_i\sim \text{Exp(2)}$ and summing $\Delta t_i$ to get event times, let's sample the inter event distances in the cumulative space as $\Delta y_i \sim \text{Exp(1)}$.

Then calculate $y_i = \Sigma_{i' = 1}^{i} \Delta y_{i'}$

Each $y_i$ can be inverted with $d(t)$ to get $t_i = d(y_i)$. 

Then the inter event times can be calculated as $\Delta t_i = t_i - t_{i-1}$, where $t_0 = 0$ for the start time.

When an event causes the rate function to change after some time $`t^{\text{*}}`$, rather than rejecting events after $`t^{\text{*}}`$ with probability 0.3, we can instead re-evaluate what each value $y_i$ in the cumulative space maps to in the time space for events $t_i$ scheduled to happen after $`t^{\text{*}}`$. Since the rate function has changed after $`t^{\text{*}}`$, for those events we calculate $t_i$ with the inverse of the new cumulative rate function, $`t_i = d^{\text{*}}(y_i)`$, and calculate the inter event times as $\Delta t_i = t_{i} - t_{i-1}$. Reschedule events with the new event times. Note that our inverse function $`d^{\text{*}}(t)`$ will be equal to $d(t)$ before $`t^{\text{*}}`$.

Let's look back at our face mask wearing scenario and walk through what this means for forecasting. Imagine that at some time $t^{\text{*}} \in [0, 1]$ an infectious agent begins wearing an imperfect mask that is effective at preventing 30% of infections. Like the previous examples the rate function $r(t)$ is 2 absent any interventions. Now we write

```math
r(t) = \begin{cases}
2, & t \in [0, t^{\text{*}}] \\
1.4, & t \in [t^{\text{*}}, 1] \\
0, & \text{else}
\end{cases}
```

then the cumulative rate function is

```math
c(t) = \begin{cases}
2t, & t \in [0, t^{\text{*}}] \\
2t^{\text{*}} + 1.4(t-t^{\text{*}}) = 1.4t+0.6t^*,  & t \in [t^{\text{*}}, 1]\\
1.4 + 0.6t^{\text{*}}, & \text{else}
\end{cases}
```

and the inverse function will be

```math
d(y) = \begin{cases}
y/2, & y \in [0, 2t^{\text{*}}] \\
\frac{y-2t^{\text{*}}}{1.4}+t^{\text{*}} = \frac{y - 0.6t^{\text{*}}}{1.4}, & y \in [2t^{\text{*}}, 1.4+0.6t^{\text{*}}] \\
\text{undefined}, & \text{else}
\end{cases}
```

Before $t^{\text{*}}$ would have been determined, we have done the following:

1. Iteratively sample $\Delta y_i \sim \text{Exp(1)}$
2. Calculate $y_i$ as $y_i = \Sigma_{i' = 1}^{i} \Delta y_{i'}$
3. Inverted $y_i$ with $d(y_i)$ to get $t_i$. Since $d(y_i) = t_i$. Once $y_i > 2$, stop sampling $\Delta y_i$.
4.  The inter event times are $\Delta t_i = t_i - t_{i-1}$.

Once we know what $`t^{\text{*}}`$ is, for every value $`t_i > t^{\text{*}}`$, do the following:

5. Map $t_i$ back to $y_i$ using the original cumulative function $c(t) = 2t$. In this example, $c(t) = 2t$, thus $y_i = 2t_i$.
6. Invert $y_i$ with $`d(y) = \frac{t - 0.6t^{\text{*}}}{1.4}`$ to get the new forecasted times for infection attempts, i.e., $`t_i^{\text{new}} = d^{\text{*}}(y_i) = \frac{y_i - 0.6t^{\text{*}}}{1.4}`$.
7. The forecasted times are now $`t_1, ..., t_i^{\text{new}}`$ for all $t \in [0, 2]$. Once $`t^{\text{new}}_i > 2`$, stop inverting $y_i$ to get $`t^{\text{new}}_i`$ since the infectious agent is no longer infectious when $t > 2$. 
8.  The inter event times are $\Delta t_i = t_i - t_{i-1}$, switching to $`\Delta t^{\text{new}}_i = t^{\text{new}}_i - t^{\text{new}}_{i-1}`$ when $`t_i > t^{\text{*}}`$. More specifically, the first inter event time calculated using the new inverse cumulative will be $`\Delta t^{\text{new}}_i = t^{\text{new}}_i - t_{i-1}`$, where $`t^{\text{new}}_i \geq t^{\text{*}}`$ and $`t_i < t^{\text{*}}`$.

For example, suppose we had forecasted $`\Delta y_i = \{0.11, 0.23, 0.74, 0.38, 0.25, 0.19, 0.53\}`$ by sampling $\Delta y_i \sim \text{Exp(1)}$. Then we would have:

1. $`\Delta y_i = \{ 0.11, 0.23, 0.74, 0.38, 0.25, 0.19, 0.53\}`$
2. $`y_i = \{ 0.11, 0.34, 1.08, 1.46, 1.71, 1.90 , 2.43\}`$
3. Invert $y_i$ with $d(y) = y/2$ to get $t_i$. In this trivial example, $`t_i = \{ 0.055, 0.17, 0.54, 0.73, 0.855, 0.95, 1.215\}`$. The last event time is past $t = 1$, so we drop it from the forecasted infection attempt times. 
4. The inter event times are $`\Delta t_i = \{ 0.055, 0.115, 0.37, 0.19, 0.125, 0.095 \}`$. 
5. Imagine that before $t = 0.75$, another random event lets us know that $t^{\text{*}} = 0.75$. Instead of accepting $t_5 = 0.855$ and $t_6 = 0.95$, we go back to the cumulative space for $y_5 = 1.71$ and $y_6 = 1.90$. 
6. Since an intervention has changed the rate and cumulative rate functions past $t = 0.75$, we'll use the modified cumulative rate function to reassess the forecasted times. Now our new $t_5$ is $t_5 = \frac{y_5 - 0.6 * 0.75}{1.4} = \frac{1.26}{1.4} = 0.9$. However, $y_6>1.4+0.6t^{\text{*}}$, so there is no $t_6$ (or equivalently $t_6 > 1$, past which our infectious agent is no longer infectious).
7. Our forecasted times are now $`t_i = \{ 0.055, 0.17, 0.54, 0.73, 0.9\}`$
8. Our inter event times are now $`\Delta t_i = \{ 0.055, 0.115, 0.37, 0.19, 0.17\}`$

Conceptually, we have found the amount of cumulative "shedding" at which a transmission occurs.  Until the mask is put on, there is no effect.  Once the mask is put on, the cumulative shedding decreases, and so we must find the times at which the cumulative shedding hits the pre-calculated thresholds.  This delays $t_5$, and the cumulative shedding never reaches the threshold that would give $t_6$.

*Observation:*
If the person took an antiviral that changes their new infectiousness, you could assume this is the person’s new intrinsic infectiousness. Sometimes, you may *have* to:
Imagine viral rebound. We have to reschedule if it *increases* future transmissibility. In this case we might expect that forecasted events happen earlier than they were originally forecasted. The approach of time scaling and forecasting events through the cumulative rate function is a straight forward approach to rescheduling events in this scenario.

## Distributing Infection Hazard

Using the example above of a person who is uniformly infectious on $`[0, 1]`$,
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

1. Forecast?
   1. Intrinsic is $`\leq 2`$ → At Home infectiousness is 4. At Work infectiousness is 8.
   2. A forecast could just be 8.
   3. Alternatively, a forecast could be 4 for time 0..1 and 8 for time 1..2.

2. Assume we use option (b). We draw a time $`t \sim \text{Exp}(1)/8`$.
   1. What if we used (c)? Cumulative hazard: $`4t`$ for 0..1; $`4 + 8t`$ for 1..2. Inverse: $`\frac{x}{4}`$ for $`x \in [0, 4]`$; $`\frac{x - 4}{8}`$ for $`x \in [4, 8]`$. Infinite otherwise.

3. Assume we get $`t = 0.5`$ from the above. Make a plan for an attempt at 0.5, with forecast 8.
4. At time 0.5: Current infectiousness?
   1. If at home: How many people are at home? Say it’s only 1 other. Current infectiousness is only 1. Accept with probability $`\frac{1}{8}`$.
   2. If at work: How many people are at work? Say it’s only 3 other. Current infectiousness is only 3. Accept with probability $`\frac{3}{8}`$.

5. Infect if contact is susceptible.

### Q: What if people can wear a facemask?

---

### Negative Binomial Branching Process Model

Q: How would we simulate a *negative binomial* branching process model with an arbitrary but fixed generation interval and density $`f(t)`$?

Q: What is the Intrinsic infectiousness?

A: **Claim:** For each person, it is proportional to $`f(t)`$.

Q: Suppose everyone had intrinsic infectiousness $`\alpha \cdot f(t)`$. What is the offspring distribution?

A: Poisson with mean $`\alpha`$!

Q: How to make it negative binomial?

A: It is a Poisson-Gamma mixture! So, choose $`\alpha`$ for each person at infection time according to a Gamma distribution.
