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

## Multiple Contacts

Imagine a scenario where you have an infectious individual with $`n`$ contacts. For each contact, we assume that infection events occur according to the same time varying hazard rate $`\lambda(t)`$. 

We want to distinguish two quantities:

1. the $`\textbf{total number of infection attempts}`$ by time $`t`$, and
2. the $`\textbf{number of distinct contacts infected}`$ by time $`t`$

Recall that these are different because the same contact can experience multiple infection attempts but can only become infected once. 


### Start with one contact

Let $`N_i(t)`$ be the number of infection attempts experienced by contact $`i`$ between time $`0`$ and time $`t`$. 

Assume these events follow a Poisson process with instantaneous hazard rate $`\lambda(t)`$. The number of events by time $`t`$
 is Poisson distributed with mean equal to the cumulative rate. Therefore, 

```math
N_i(t) \sim \operatorname{Poisson}(\Lambda(t))
```

or equivalently, 

```math
N_i(t) \sim \operatorname{Poisson}\left(- \int_0^t \lambda(s)\, ds \right)
```


Thus the probability of any particular number of infection attempts is 

```math
P(N_i(t) = k) = \frac{e^{-\Lambda(t)} \Lambda(t)^k}{k!}, k = 0, 1, 2, ...
```

and 

```math
\mathbb{E}[N_i(t)] = \Lambda(t)
```

So $`\Lambda(t)`$ is the expected number of infection attempts experienced by one contact by time $`t`$. 

### Probability that one contact is infected

Suppose the first infection event is sufficient to infect the contact. The contact remains uninfected by time $`t`$ exactly when it has experienced zero infection attempt events. 

```math
P(\text{not infected by time } t) = P(N_i(t) = 0)
```

Since $`N_i(t)`$ is Poisson distributed, then using the expression above we have, 

```math
P(N_i(t) = 0) = \frac{e^{-\Lambda(t)} \Lambda(t)^0}{0!} = e^{-\Lambda(t)}
```

Since a contact can either be infected or not by time $`t`$, we can say that the probability of infection by time $`t`$ is

```math
P(\text{infected by time }t) = 1 - e^{-\Lambda(t)} = p_t
```

Which is the same relationship from our definition of the hazard function and the survival function above. 

### Infection status is a Bernoulli random variable

For a contact $`i`$ we define

```math
X_i(t) =
\begin{cases}
    1, \text{ if contact is infected by time } t, \\
    0, \text{ otherwise}
\end{cases}
```

We just established that 

```math
P(X_i(t) = 1) = p_t
```

and therefore

```math
P(X_i(t) = 0) = 1 - p_t
```

This is exactly the definition of a Bernoulli random variable, so

```math
X_i(t) \sim \operatorname{Bernoulli}(p_t)
```

### Extending to $`n`$ contacts

Now imagine a scenario where an infected person has $`n`$ contacts they could infect and each infection process is independent. We therefore have

```math
N_1(t), N_2(t), ..., N_n(t)
```

where $`N_i(t) \sim \operatorname{Poisson}(\Lambda(t))`$ .

Now we can derive the population-level quantities. 

Let 

```math
N(t) = \sum_{i=1}^n N_i(t)
``` 

be the total number of infection attempts across all contacts by time $`t`$. 

Independent Poisson random variables add linearly so we can say that

```math
N(t) \sim \operatorname{Poisson}\left( \sum_{i=1}^n \Lambda(t) \right)
```

Since every contact has the same cumulative hazard,

```math
\sum_{i=1}^n \Lambda(t) = n\Lambda(t)
```

Therefore 

```math
N(t) \sim \operatorname{Poisson}(n\Lambda(t))
```

or we can also write this as 

```math
N(t) \sim \operatorname{Poisson}\left(n \int_0^t \lambda(s)\, ds\right)
```

### Number of distinct contacts infected

Let

```math
X(t) = \sum_{i=1}^n X_i(t)
```

be the number of distinct contacts infected by time $`t`$. 

For each contact, $`X_i(t) \sim \operatorname{Bernoulli}(p_t)`$. Again, assuming the infection processes of the contacts are independent, then $`X_1(t), X_2(t), ..., X_n(t)`$ are independent Bernoulli random variables with the same probability $`p_t`$. The sum of $`n`$ independent Bernoulli random variables with the same probability is Binomial, thus $`X(t) \sim \operatorname{Bernoulli}(p_t)`$. 

Why?

Let's consider one particular outcome in which exactly $`k`$ contacts are infected by time $`t`$. The probability of one way to get this outcome with $`n`$ contacts is

```math
p_t^k (1 - p_t)^{n-k}
```

i.e., $`k`$ contacts get infected and $`n-k`$ contacts do not get infected by time $`t`$, and since infection of each contact is independent we can multiply the probability of each contact's state.

But there are many ways for $`k`$ contacts to be infected. We can choose the $`k`$ infections to happen to any of the $`n`$ contacts. The number of way to choose which contacts are infected is

```math
n\choose k
```

So the probability of infecting any $`k`$ contacts can be defined as 

```math
P(X(t) = k) = {n\choose k} p_t^k (1 - p_t)^{n-k}
```

This is exactly the probability mass function of a Binomial random variable. 

Therefore, for $`n`$ contacts we have

```math
N(t) \sim \operatorname{Poisson}\left( n\int_0^t \lambda(s) \, ds \right)
```
```math
X(t) \sim \operatorname{Binomial}\left(n, 1 - \operatorname{exp}\left( - \int_0^t \lambda(s) \, ds \right) \right)
```

to describe the number of infection attempts and the number of infections caused by an infectious person up to time $`t`$. 

For person 1, let's look again at the probability that they become infected when the infectious person has $`n`$ contacts. The infectious person will generate a total number of infection attempts $`k`$ by time $`t`$ according to a Poisson distribution with mean $`n\Lambda(t)`$, 

```math
P(N(t) = k) = \frac{(n\Lambda(t))^ke^{-n\Lambda(t)}}{k!}
```

Each attempt targets person 1 with probability $`\frac{1}{n}`$, so person escapes a given attempt with probability $`1 - \frac{1}{n} = \frac{n-1}{n}`$. Given $`k`$ attempts, the probability of person 1 not being infected by time $`t`$ is $`\left(\frac{n-1}{n}\right)^k`$. 

Summing over all possible values of $`k`$ we have,

```math
P(\text{person 1 not infected}) = \sum_{k=0}^\infty \left( \frac{n-1}{n}\right)^k \frac{ \left(n\Lambda(t) \right)^k e^{-n\Lambda(t)} }{k!}
```

Rearranging we get

```math
\begin{align*}
P(\text{person 1 not infected}) 
&= e^{-n\Lambda(t)} \sum_{k = 0}^\infty \left( \frac{n-1}{n}\right)^k \frac{\left( n\Lambda(t)\right)^k}{k!} \\
&= e^{-n\Lambda(t)} \sum_{k = 0}^\infty \left( n - 1\right)^k \frac{n^k}{n^k} \frac{\left(\Lambda(t)\right)^k}{k!} \\
&= e^{-n\Lambda(t)} \sum_{k = 0}^\infty \frac{\left( (n - 1)\Lambda(t) \right)^k}{k!}
\end{align*}
```

Recall the Taylor series expansion of $`e^x`$:

```math
e^x = \sum_{n = 0}^\infty \frac{x^n}{n!}
```

Therefore we have


```math
\begin{align*}
P(\text{person 1 not infected}) 
&= e^{-n\Lambda(t)} e^{(n-1)\Lambda(t)} \\
&= e^{-\Lambda(t)}
&= 1 - p_t
\end{align*}
```

Then the probability of person 1 being infected is $`p_t`$, just as we found before. 


<!-- 
Suppose infection events for a person follow a Poisson process with rate function $`\lambda(t)`$. This means that you are scheduling a Poisson process with rate function $`\lambda(t)`$. 

Let $`N(t)`$ be the number of events by time $`t`$. Then -->
<!-- 
```math
N(t) \sim \operatorname{Poisson}\left( \int_0^t \lambda(s)\, ds \right)
``` -->
<!-- 
The person is infected by time $`t`$ if at least one event has occurred:

```math
P(\text{infected by time } t) = P(N(t) \geq 1)
```
For a Poisson random variable, 

```math
P(N(t) = 0) = \operatorname{exp}\left ( -\int_0^t \lambda(s)\, ds \right)
``` -->

<!-- Therefore, 

```math
\begin{align*}
P(N(t) \geq 1) &= 1 - P(N(t) = 0) \\
 &= 1 - \operatorname{exp}\left(-\int_0^t \lambda(s)\, ds \right) \\
 &= p_t
\end{align*}
```

So at each event you attempt to infect a given person and by time $`t`$ the person will be infected with probability $`p_t`$. -->

<!-- ## Proposition

Imagine you have a setting where each person has $`n`$ contacts and there is an infectious person in the setting. We assume that each contact independently experiences the same transmission hazard $`\lambda(t)`$. Infection attempts in this scenario can be modeled as a Poisson process with rate function $`n\lambda(t)`$.  -->

<!-- you schedule a Poisson process of infection attempts with rate function $`n \lambda(t)`$. -->

<!-- Then the number of people infected by time $`t`$ is given by $`\text{Binomial}(n, p_t)`$.

## Working this out for a single contact in this setting

Let's work this out for a single contact in this setting, call them person 1.

For each contact $`i`$, we can define

```math
X_i(t) = 
\begin{cases}
    1, \text{ if contact } i \text{ is infected by time } t, \\
    0, \text{ otherwise}
\end{cases}
``` -->

<!-- From the one-person result above, we know that 

```math
P(X_1(t) = 1) = p_t
```

Since person 1 can either be infected or not after infection attempts, 

```math
P(X_1(t) = 0) = 1 - p_t
``` -->

<!-- This is exactly the definition of a Bernoulli random variable. A random variable $`X`$ has a Bernoulli distribution with parameter $`p`$ when

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
``` -->

<!-- 
Because the infection attempts against each contact are independent, it follows that $`X_i(t)`$ for each contact is also an independent Bernoulli random variable with the same transmission probability $`p_t`$. Therefore their sum, 

```math
X(t) = \sum_{i=1}^n X_i(t)
```

is the sum of $`n`$ independent Bernoulli random variables with the same probability.  -->

<!-- For $`n`$ contacts, consider the particular outcome in which exactly $`k`$ contacts are infected. The probability of one way to get this outcome is

```math
p_t^k (1 - p_t)^{n-k}
```

But there are multiple ways for $`k`$ contacts to be infected -  we can choose the $`k`$ infections to happen to any of the $`n`$ contacts. The number of ways to choose which $`k`$ contacts are infected out of all $`n`$ contacts is 

```math
n\choose k
``` -->

<!-- Each of these possible ways has the same probability

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
``` -->


<!--  
Let's say $`c = \int_0^t \lambda(s)\, ds`$, which is the cumulative hazard for given a single contact.

For setting with `n` contacts, the number of attempts by time $`t`$ will be $`\text{Poisson}(\text{rate} = nc)`$. -->
<!-- 
So if $`k`$ is the number of attempts, then that has Poisson probability

```math
\frac{(nc)^k e^{-nc}}{k!}
```

Each attempt has probability $`\frac{1}{n}`$ of being directed at person 1 (they escape being infected with probability $1-\frac{1}{n}$). The probability that person 1 has *not* been infected given there are $`k`$ attempts is $`\left(\frac{n-1}{n}\right)^k`$.

So the total probability that person 1 has not been infected is (i.e., summing over attempts to infinity)

```math
\sum_k \left(\frac{n-1}{n}\right)^k \frac{(nc)^k e^{-nc}}{k!}
``` -->

<!-- By rearranging:

```math
= \sum_k \frac{\big((n-1)c\big)^k e^{-nc}}{k!}
```

but remember that $`e^x = \sum_k \frac{x^k}{k!}`$, so the above equals

```math
e^{(n-1)c} e^{-nc}
```

```math
= e^{nc - c - nc} = e^{-c} = 1 - p_t
``` -->
