# 2.0 Hazard with Multiple Contacts

## Recall

If $`\lambda(t)`$ is the hazard of transmission, then $`p_t`$, the probability of infection by time $`t`$, is

```math
p_t = 1 - \exp\left(-\int_0^t \lambda(s)\, ds\right)
```

## Fact

If you schedule a Poisson process with rate function $`\lambda(t)`$, and at each event you attempt to infect a given person, then the person will be infected with probability $`p_t`$.

## Another fact

Imagine you have a setting with $`n`$ people and you schedule a Poisson process of infection attempts with rate function $`n \lambda(t)`$.

Then the number of people infected by time $`t`$ is actually $`\text{Binomial}(n, p_t)`$.

## Working this out for a single person in this setting

Let's work this out for a single person in this setting, say person 1.

Let's say $`c = \int_0^t \lambda(s)\, ds`$.

The number of attempts by time $`t`$ will be $`\text{Poisson}(\text{rate} = nc)`$.

So if $`k`$ is the number of attempts, then that has Poisson probability

```math
\frac{(nc)^k e^{-nc}}{k!}
```

The probability that person 1 has *NOT* been infected given there are $`k`$ attempts is $`\left(\frac{n-1}{n}\right)^k`$.

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
