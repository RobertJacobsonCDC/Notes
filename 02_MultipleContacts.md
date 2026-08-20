# 2.0 Hazard with Multiple Contacts

## Recall

Let $`\lambda(t)`$ be the instantaneous hazard rate of transmission at time $`t`$ for one susceptible person. Define the cumulative hazard rate function through time $t$ as

```math
\Lambda(t) = \int_0^t \lambda(s)\, ds
```

The probability that the person has yet to be infected by time $`t`$ is 

```math
S(t) = e^{-\Lambda(t)} = \exp\left( \int_0^t \lambda(s)\,  ds\right)
```

Then $`p_t`$, the probability of infection by time $`t`$, is

```math
p_t = 1 - \exp\left(-\int_0^t \lambda(s)\, ds\right)
```

This is the basic relationship between the transmission hazard and the probability of infection.

## Lemma

Suppose infection events for a person follow a Poisson process with rate function $`\lambda(t)`$. This means that you are scheduling a Poisson process with rate function $`\lambda(t)`$. 

Let $`N(t)`$ be the number of events by time $`t`$. Then

```math
N(t) \sim \operatorname{Poisson}\left( \int_0^t \lambda(s)\, ds \right)
```

The person is infected by time $`t`$ if at least one event has occurred:

```math
P(\text{infected by time } t) = P(N(t) \geq 1)
```
For a Poisson random variable, 

```math
P(N(t) = 0) = \operatorname{exp}\left ( -\int_0^t \lambda(s)\, ds \right)
```

Therefore, 

```math
\begin{align*}
P(N(t) \geq 1) &= 1 - P(N(t) = 0) \\
 &= 1 - \operatorname{exp}\left(-\int_0^t \lambda(s)\, ds \right) \\
 &= p_t
\end{align*}
```

So at each event you attempt to infect a given person and by time $`t`$ the person will be infected with probability $`p_t`$.

## Proposition

Imagine you have a setting where each person has $`n`$ contacts and there is an infectious person in the setting. We assume that each contact independently experiences the same transmission hazard $`\lambda(t)`$. Infection attempts in this scenario can be modeled as a Poisson process with rate function $`n\lambda(t)`$. 

<!-- you schedule a Poisson process of infection attempts with rate function $`n \lambda(t)`$. -->

Then the number of people infected by time $`t`$ is given by $`\text{Binomial}(n, p_t)`$.

## Working this out for a single contact in this setting

Let's work this out for a single contact in this setting, call them person 1.

For each contact $`i`$, we can define

```math
X_i(t) = 
\begin{cases}
    1, \text{ if contact } i \text{ is infected by time } t, \\
    0, \text{ otherwise}
\end{cases}
```

From the one-person result above, we know that 

```math
P(X_1(t) = 1) = p_t
```

Since person 1 can either be infected or not after infection attempts, 

```math
P(X_1(t) = 0) = 1 - p_t
```

This is exactly the definition of a Bernoulli random variable. A random variable $`X`$ has a Bernoulli distribution with parameter $`p`$ when

```math
P(X = 1) = p
```

and

```math
P(X = 0) = 1 - p
```

Thus, we can say for person 1,

```math
X(t) \sim \operatorname{Bernoulli}(p_t)
```


Because the infection attempts against each contact are independent, it follows that $`X_i(t)`$ for each contact is also an independent Bernoulli random variable with the same transmission probability $`p_t`$. Therefore their sum, 

```math
X(t) = \sum_{i=1}^n X_i(t)
```

is the sum of $`n`$ independent Bernoulli random variables with the same probability. 

For $`n`$ contacts, consider the particular outcome in which exactly $`k`$ contacts are infected. The probability of one way to get this outcome is

```math
p_t^k (1 - p_t)^{n-k}
```

But there are multiple ways for $`k`$ contacts to be infected -  we can choose the $`k`$ infections to happen to any of the $`n`$ contacts. The number of ways to choose which $`k`$ contacts are infected out of all $`n`$ contacts is 

```math
n\choose k
```

Each of these possible ways has the same probability

```math
p_t^t (1 - p_t)^{n-k}
```

So the probability of infecting $`k`$ contacts can be defined as

```math
P(X(t) = k) = {n\choose k} p_t^k (1 - p_t)^{n-k}
```

This is exactly the probability mass function of a Binomial random variable. Therefore, 

```math
X(t) \sim \operatorname{Binomial}(n, p_t).
```


 
Let's say $`c = \int_0^t \lambda(s)\, ds`$, which is the cumulative hazard for given a single contact.

For setting with `n` contacts, the number of attempts by time $`t`$ will be $`\text{Poisson}(\text{rate} = nc)`$.

So if $`k`$ is the number of attempts, then that has Poisson probability

```math
\frac{(nc)^k e^{-nc}}{k!}
```

Each attempt has probability $`\frac{1}{n}`$ of being directed at person 1 (they escape being infected with probability $1-\frac{1}{n}$). The probability that person 1 has *not* been infected given there are $`k`$ attempts is $`\left(\frac{n-1}{n}\right)^k`$.

So the total probability that person 1 has not been infected is (i.e., summing over attempts to infinity)

```math
\sum_k \left(\frac{n-1}{n}\right)^k \frac{(nc)^k e^{-nc}}{k!}
```

By rearranging:

```math
= \sum_k \frac{\big((n-1)c\big)^k e^{-nc}}{k!}
```

but remember that $`e^x = \sum_k \frac{x^k}{k!}`$, so the above equals

```math
e^{(n-1)c} e^{-nc}
```

```math
= e^{nc - c - nc} = e^{-c} = 1 - p_t
```
