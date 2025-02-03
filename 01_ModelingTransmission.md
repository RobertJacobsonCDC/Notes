# 1.0 Parameterizing and Modeling Transmission

## Simple Poisson Process

Assume a person has intrinsic infectiousness of 1 person per day for time $[0, 2]$.

By definition, I mean:
$\displaystyle \frac{d}{dt} P(\text{Transmit at time } t \mid \text{Have not transmitted by time } t)  = 1 $.
This means:
$$\displaystyle \lim_{\varepsilon \to 0} \frac{P(\text{transmit in } [t, t + \varepsilon] \mid \text{have not transmitted by time } t)}{\varepsilon}. $$


### Q: How would you simulate the probability and time of infection for a single susceptible contact?

**Answer:** Sample $x \sim \text{Exp}(1)$. If $x \lt  2$, then the person is infected at time $x$. Otherwise, the person is not infected.

### Q: How would you simulate a Poisson process with this same rate function?

$$
r(t) = \begin{cases}
1, & t \in [0, 2] \\
0, & \text{else}
\end{cases}
$$

**Option 1:** First compute the number of events, then distribute in time.
Sample from a $\text{Poisson}(2) \to n$. Then sample $n$ times uniformly in $[0, 2]$.

*Intuition:* Case $n = 1$: We are assuming $x_1 \lt  2$ AND $x_1 + x_2 > 2$.

**Option 2:** Iteratively sample $x_i \sim \text{Exp}(1)$. Then the times are $x_1, x_1 + x_2, x_1 + x_2 + x_3,\dots$ up to when they exceed 2.

### Observation:
If we simulate from a Poisson process of rate $r(t)$, then the probability of having one event and the time of that event matches the intrinsic infectious process.

## Rejection Sampling and Time Scaling

Consider how to handle changes:
Suppose at a random time $t \in [0, 2]$, this person will wear a perfect facemask.

### Q: How does this change intrinsic transmissibility, and how do we simulate?

**Option 1:** Rejection sample. You sample according to Option 2. Then at the time of the purported event, if it is after the facemask time $t$, you reject it (ignore it).

What if the facemask isn’t perfect? Say it only prevents transmission 30% of the time.

**Answer:** You still schedule it normally, and then reject with probability 0.3.

*Observation:*
This 30% is not the reduction in the probability of infecting a single contact. It is the reduction in the probability of a single infection attempt succeeding.

**Option 2:** When facemask time comes, you cancel the existing plan and schedule based on the new rate (of 0.7 up to time 2).
The answer is to sample from $\text{Exp}(1) / 0.7$, which is equivalent to an Exponential distribution of rate 0.7 — this is referred to as "Time rescaling."

*Observation:*
If the person took an antiviral that changes their new infectiousness, you could assume this is the person’s new intrinsic infectiousness. Sometimes, you may *have* to:
Imagine viral rebound. We have to reschedule if it *increases* future transmissibility.

## Distributing Infection Hazard

Using the example above of a person who is uniformly infectious on $[0, 2]$,
suppose an infected person lives in a house with 2 other people, and they stay in that house for all time during the course of the infection.

### How should we model the infection process of those two people in terms of an inhomogeneous Poisson process?

**Possibility 1:** The hazard of infecting each person is equal and half of the intrinsic infectiousness. The rate of the total infection process is, in this case, just the intrinsic infectiousness (halved and then added for each person).

**Possibility 2 :** The hazard of infecting each person is equal to the intrinsic infectiousness. The rate process of the total infection process is twice the intrinsic infectiousness.

In the real world we observe something between these two extremes.

### Example from History:
We assumed there was a magic parameter $\alpha$ and that in a setting with $n$ people (total, including the source of infection), the rate of the total infectiousness process was:

$$
\text{Intrinsic Infectiousness} \times (n - 1)^\alpha \quad \text{where} \quad 0 \leq \alpha \leq 1,
$$

which interpolates between the two options.


## Modeling Transmission with Time-Varying Factors

Suppose this person, who is in a house with 2 people, is going to wear a facemask at random times throughout the day.

### Claim:
Let $r(t)$ be the intrinsic infectiousness for this person in the absence of any interventions. You schedule infections according to $2^\alpha \cdot r(t)$ and reject them with probability 0.3 if the person is wearing a facemask at that time and otherwise accept. This works because $r(t) \mapsto 2^\alpha \cdot r(t)$ is linear in $r(t)$.

However, total infectiousness in this setting could be some function $f_2(r(t))$ for some arbitrary monotonic positive function. In this case, you would accept with probability:

$$ \frac{f_2(0.7 \cdot r(t))}{f_2(r(t))} $$

### Linearity is Convenient:
Total infectiousness can be modeled as:

$$ \text{Total infectiousness} = \left( \text{fraction of time in home} \cdot \text{home multiplier} + \text{fraction of time in work} \cdot \text{work multiplier} \right) \cdot r(t) $$

Suppose a person spends half of their time at home and half at work. Home has 3 people, work has 9 people. Assume $\alpha = 1$ as above.
The total infectiousness would be
$\left( 0.5 \cdot 2 + 0.5 \cdot 8 \right) \cdot r(t)$.

### How Would We Handle Telework?
From the perspective of the source, you adjust the fraction of time.

**Proposal:** Schedule infection attempts at a rate of
$(0.5 \cdot 2 + 0.5 \cdot 8) \cdot r(t)$
and record the expected total infectiousness at that time in the plan. Call that $R$.
At the time of the attempt, compute the person’s total infectiousness in each setting. Let $S$ be the sum of these infectiousnesses.
You accept the proposed attempt with probability
$S/R$.
Then, you choose the setting proportional to their current total infectiousness. After that, you pick the contact in that setting.

Assume a person is either at home or at work at any point in time, but you don’t know the schedule.
At home, their total infectiousness would be $2 \cdot r(t)$. At work, it would be $8 \cdot r(t)$.

**Proposal:** Schedule based on $8 \cdot r(t)$. Proceed as above.

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

$$
\text{intrinsic infectiousness} =
\begin{cases}
1, & 0 \leq t \lt  1 \\
2, & 1 \leq t \lt  2 \\
0, & t \geq 2
\end{cases}
$$

For a specific person who was just infected, assume:
- Their Home setting has up to 3 people in it (including themselves).
- Their Work setting has up to 5 people in it (including themselves).

Assume infection risk applies independently to people in all settings ($\alpha = 1$ as above).

**Loop:**

1. Forecast?
   1. Intrinsic is $\leq 2$ → At Home infectiousness is 4. At Work infectiousness is 8.
   2. A forecast could just be 8.
   3. Alternatively, a forecast could be 4 for time 0..1 and 8 for time 1..2.

2. Assume we use option (b). We draw a time $t \sim \text{Exp}(1)/8$.
   1. What if we used (c)? Cumulative hazard: $4t$ for 0..1; $4 + 8t$ for 1..2. Inverse: $\frac{x}{4}$ for $x \in [0, 4]$; $\frac{x - 4}{8}$ for $x \in [4, 8]$. Infinite otherwise.

3. Assume we get $t = 0.5$ from the above. Make a plan for an attempt at 0.5, with forecast 8.
4. At time 0.5: Current infectiousness?
   1. If at home: How many people are at home? Say it’s only 1 other. Current infectiousness is only 1. Accept with probability $\frac{1}{8}$.
   2. If at work: How many people are at work? Say it’s only 3 other. Current infectiousness is only 3. Accept with probability $\frac{3}{8}$.

5. Infect if contact is susceptible.

### Q: What if people can wear a facemask?

---

### Negative Binomial Branching Process Model

Q: How would we simulate a *negative binomial* branching process model with an arbitrary but fixed generation interval and density $f(t)$?

Q: What is the Intrinsic infectiousness?

A: **Claim:** For each person, it is proportional to $f(t)$.

Q: Suppose everyone had intrinsic infectiousness $\alpha \cdot f(t)$. What is the offspring distribution?

A: Poisson with mean $\alpha$!

Q: How to make it negative binomial?

A: It is a Poisson-Gamma mixture! So, choose $\alpha$ for each person at infection time according to a Gamma distribution.
