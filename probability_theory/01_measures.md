# Probability Theory is a Subfield of Measure Theory

We frequently see definitions like the following:

> ```math
> P(a \lt X\lt b):=\int_a^bf_X(x) dx
> ```

and

> **Cumulative Density Function (CDF)**: A function giving the probability that a continuous random variable is less than or equal to a specific value, $F_X(x) := P(X \leq x) = \int_{-\infty}^x f_X(w) dw$.

While these definitions are useful in their own right, there is a deeper reason why they are given in terms of intervals: they define a measure. As we will see below:

> The definition of a measure on intervals (or equivalently, half intervals) _completely determines_ what the measure is on _every_ measurable set.

## Measures And The Borel Measure on $\mathbb{R}$

A measure is a way to define a notion of the "size" of a subset of some ambient space. The classic example is the _Borel measure_ $B$, which gives us a way of measuring a very large family of subsets of $\mathbb{R}$, called _Borel sets_ (or _Borel measurable sets_), which we denote $\mathcal{B}$. The Borel measure is exactly the measure that corresponds to _length_ when applied to open intervals: $B((a, b)) := b-a$ for $a, b \in \mathbb{R}$. If we want $B$ to be a useful and consistent notion of "size," then we obviously want:

1. $B(\varnothing) = 0$.
2. $B(A) \geq 0$ for any $A\in \mathcal{B}$ (for any $A$ that can be measured).
3. For any sequence $\{A_n\}_{n=1}^\infty$ of pairwise disjoint measurable sets, ${\displaystyle B {\left(\bigcup _{n=1}^{\infty }A_{n}\right)}=\sum _{n=1}^{\infty }B (A_{n})}$.

It turns out these are all the axioms we need for any measure. We can then "construct" the entire family of measurable sets by starting with some "generating set"—the open intervals in the case of the Borel sets—and then doing either one of the following:

1. Construct the family of measurable sets "intrinsicly" (from the inside out) by iteratively taking complements and countable unions, at each stage extending the measure itself uniquely to the new sets so constructed using the measure axioms (1)-(3) above.
2. Construct the family of measurable sets "extrinsically" (from the outside in) by taking the intersection of all $\sigma$-algebras that contain the generating set, and show that the result is the unique smallest $\sigma$-algebra containing the generating set.

It turns out that either way you get a $\sigma$-algebra of sets, which by definition is a family of sets that:
1. contains the original ambient set (in the Borel case, $\mathbb{R}\in\mathcal{B}$);
2. is closed under set complement;
3. is closed under countable unions.

Notice that how the measure is defined on the generating set _completely determines_ the measure on every measurable set.

## A Probability Space Is A Measure Space

A measure space $(\Omega, \mathcal{F}, \mu)$ is defined by three components:

1. A nonempty set $\Omega$.
2. A $\sigma$-algebra $\mathcal{F}$ of subsets of $\Omega$, the measurable sets.
3. A measure $\mu \colon \mathcal{F} \to \mathbb{R}\cup\{+\infty\}$ satisfying the three measure axioms given in the previous section.

A probability space $(\Omega, \mathcal{F}, P)$ is just a measure space ($P=\mu$) for which $P(\Omega) = 1$. 

## Why all this technical machinery?

Why all this technical machinery? Because it turns out we can't meaningfully define a reasonable notion of "length" / "area" for _arbitrary_ subsets.

> **Banach–Tarski:** For any arbitrary bounded $U$ and $V$, there is a partition of $U$ into finite disjoint $\{E_n\}_{n=1}^N$ s.t. rigid rotations and translations of $\{E_n\}$ union to $V$.

Rigid motions of disjoint unions do not preserve volume. Our only hope is to restrict ourselves to some family of subsets for which we can define a meaningful measure. These formalisms do exactly that.
