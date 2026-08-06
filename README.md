# Math and Models

Notes on Mathematical Statistics and Statistical Modeling

## Contents

- [Section 0: Background Information](00_Background.md)
- [Section 1: Parameterizing and Modeling Transmission](01_ModelingTransmission.md)
- [Probability Theory](probability_theory/index.md)
  - [Probability Spaces](probability_theory/01_probability_space.md)
  - [Modeling Wait Time](probability_theory/02_wait_time.md)
  - [From Wait Time to Poisson Process](probability_theory/03_poisson_process.md)

The probability theory notes are organized as a VitePress book. The Markdown
files in the repository root remain standalone notes for now.

## Development

Install the JavaScript dependencies with:

```sh
pnpm install
```

Common VitePress tasks are available through
[`plzplz`](https://plzplz.org/):

```sh
plz docs dev       # Start the development server with hot reload
plz docs build     # Build the production site
plz docs preview   # Build and preview the production site locally
```

## Math Syntax

The notes use GitHub-compatible math delimiters. Inline expressions use dollar
signs with an inner code span, and displayed expressions use a `math` code
fence:

````markdown
Inline: $`\{E_n\}_{n=1}^N`$

```math
P_X := P \circ X^{-1}
```
````

The VitePress configuration translates both forms to MathJax during builds.
