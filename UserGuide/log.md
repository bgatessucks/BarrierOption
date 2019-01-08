# Questions

  - How to get the analytical derivatives through the multi variate `CDF` ?
  - How to define the different variations starting from one generic payoff ?

# Notation

|symbol|legend|
|------|------|
|S|underlying|
|L|lower barrier|
|U|upper barrier|
|X|strike|
|t$^{\ast}$|hitting time|

## Atom: `up-and-out; down-and-out`

$$
\text{payoff} = \max(S-X, 0)\ \mathbb{1}_{L<S<U}\ \equiv\ \text{uodo}(S, X, L, U, \{\})
$$

# Types


  - `call; up-and-out; down-and-out`: $\text{payoff}(S, U, L) = \max(S-X, 0)\ \mathbb{1}_{L<S<U}$

  - `put; up-and-out; down-and-out`: $\text{payoff}(S, U, L) = \max(X-S, 0)\ \mathbb{1}_{L<S<U}$


# The zoo

  - `double knock-in`: $S_{0}<L\ || \ U<S_{0}$, $ \max(S-X, 0)\ \mathbb{1}_{L<S(t^{\ast})<U} $    
  - `down-and-in`: $U<S_{0}$, $\text{payoff} = \max(S-X, 0)\ \mathbb{1}_{S(t^{\ast})<U}$


# Junk
