# Pieri's

This is a tool for computing the cohomology ring structure of the Grassmannian G(k,n+1) using Schubert calculus. Namely, the tool can compute the cohomology groups using Schubert classes which serve as generating cocycles in the cohomology ring. Using Pieri's formula, we can compute the product between Schubert classes.

## Documentation

For some basic examples of usage:

```python
# To define a Schubert cycle, we use the SchubertCycle class
s1 = SchubertCycle(a=[1],n=3)
# We can then multiply cycles together or add them using * and +
s1*s1
s1 + s1
```

```
σ_(1,1)+σ_2 # output of product
s1 + s1
```

We can also find the generators for the cohomology ring

```python
compute_cohomology_groups(k=2,n=3)
[[1], [], [σ_1], [], [σ_2, σ_(1,1)], [], [σ_(2,1)], [], [σ_(2,2)]]
# The output of this function is stored as a list of lists where the i-th list are the generators of the i-th cohomology class
```

While this can be read easily from the output above, if we want a prettified version of the cohomology groups with generators, we can do this easily

```python
print_cohomology_groups(k=2,n=3)
```

```
H^0 (G(2,4)) = ℤ1
H^1 (G(2,4)) = 0
H^2 (G(2,4)) = ℤσ_1
H^3 (G(2,4)) = 0
H^4 (G(2,4)) = ℤσ_2⊕ℤσ_(1,1)
H^5 (G(2,4)) = 0
H^6 (G(2,4)) = ℤσ_(2,1)
H^7 (G(2,4)) = 0
H^8 (G(2,4)) = ℤσ_(2,2)
```

Finally, we can compute the multiplication table

```python
groups = compute_cohomology_groups(k=2,n=3)
print_multiplication_table(groups)
```

```
|     *               |     1               |     σ_1             |     σ_2             |     σ_(1,1)         |     σ_(2,1)         |     σ_(2,2)         |
-----------------------------------------------------------------------------------------------------------------------------------------------------------
|     1               |     1               |     σ_1             |     σ_2             |     σ_(1,1)         |     σ_(2,1)         |     σ_(2,2)         |
|     σ_1             |     σ_1             |     σ_(1,1)+σ_2     |     σ_(2,1)         |     σ_(2,1)         |     σ_(2,2)         |     0               |
|     σ_2             |     σ_2             |     σ_(2,1)         |     σ_(2,2)         |     0               |     0               |     0               |
|     σ_(1,1)         |     σ_(1,1)         |     σ_(2,1)         |     0               |     σ_(2,2)         |     0               |     0               |
|     σ_(2,1)         |     σ_(2,1)         |     σ_(2,2)         |     0               |     0               |     0               |     0               |
|     σ_(2,2)         |     σ_(2,2)         |     0               |     0               |     0               |     0               |     0               |
```

There is also an option to output the multiplication table into a LaTeX tabular environment
```python
groups = compute_cohomology_groups(k=2,n=3)
print_multiplication_table(groups, latex=True)
```
```
\begin{tabular}{ | c | c | c | c | c | c | c | }
\hline
 $\cdot$ &  $1$ &  $\sigma_{1}$ &  $\sigma_{1,1}$ &  $\sigma_{2}$ &  $\sigma_{2,1}$ &  $\sigma_{2,2}$ \\
\hline
 $1$ &  $1$ &  $\sigma_{1}$ &  $\sigma_{1,1}$ &  $\sigma_{2}$ &  $\sigma_{2,1}$ &  $\sigma_{2,2}$ \\
\hline
 $\sigma_{1}$ &  $\sigma_{1}$ &  $\sigma_{1,1}+\sigma_{2}$ &  $\sigma_{2,1}$ &  $\sigma_{2,1}$ &  $\sigma_{2,2}$ &  $0$ \\
\hline
 $\sigma_{1,1}$ &  $\sigma_{1,1}$ &  $\sigma_{2,1}$ &  $\sigma_{2,2}$ &  $0$ &  $0$ &  $0$ \\
\hline
 $\sigma_{2}$ &  $\sigma_{2}$ &  $\sigma_{2,1}$ &  $0$ &  $\sigma_{2,2}$ &  $0$ &  $0$ \\
\hline
 $\sigma_{2,1}$ &  $\sigma_{2,1}$ &  $\sigma_{2,2}$ &  $0$ &  $0$ &  $0$ &  $0$ \\
\hline
 $\sigma_{2,2}$ &  $\sigma_{2,2}$ &  $0$ &  $0$ &  $0$ &  $0$ &  $0$ \\
\hline
\end{tabular}
```
![Output of Multiplication Table in LaTeX](./latexoutputtest.png)


The number of generators grow quickly for higher `k` and `n`, so we will not include many examples. But it should be easy to extrapolate if needed. For one, fairly useless example,

```python
s2 = SchubertCycle(a=[2],k=4,n=9)
s43 = SchubertCycle(a=[4,3],k=4,n=9)
s2*s43
```
```
σ_(4,3,2)+σ_(4,4,1)+σ_(5,3,1)+σ_(5,4)+σ_(6,3)
```

One can also subtract cycles, negate cycles, multiply cycles by integer scalars, and as alluded to before exponentiate. All of these follow the conventional pythonic style

```python
s1 = SchubertCycle(a=[1],k=2,n=3)
s1 + s1
s1 - s1
3*(s1**2)
```
```
2σ_(1)
0
3σ_(1,1)+3σ_(2)
```

## Applications

Finally, let's use the software to solve a couple enumerative geometry problems.

### Problem 1

Q: How many lines pass through the 4 distinct lines in projective 3-space, `P^3`?

To solve this, we first translate the problem to Schubert calculus. The set of lines that pass through a given line is the Schubert cycle `σ_1`. We can then take the intersection pairing on `σ_1^4` to find the number of lines passing through any four given lines (in general position).
Let's translate this to code.

```python
s1 = SchubertCycle(a=[1], k=2, n=3)
s1**4
# or
(s1*s1)*(s1*s1)
# or
s1*s1*s1*s1
```

In any case, we obtain the output

```
2σ_(2,2)
```

And the intersection pairing simply yields the scalar 2. Hence, there are 2 lines that intersect 4 arbitrary lines in `P^3`

### Problem 2
Q: How many lines pass through 6 distinct planes in projective 4-space, `P^4`?

Equivalently the problem boils down to computing `σ_1^6` in Schubert calculus since `σ_1` can be though of as the set of lines passing through a particular plane.
Let's translate this to code.

```python
s1 = SchubertCycle(a=[1],k=2,n=4)
s1**6
```

Which yields the output

```
5σ_(3,3)
```

The intersection pairing gives the scalar 5. Hence, there are 5 lines that intersect 6 arbitrary planes in `P^4`.

### Problem 3
Q: How many lines pass through 3 distinct planes and a specific point in projective 4-space, `P^4`?

Similar to problem 2, the problem boils down to computing `(σ_1^3)*(σ_3)`. We can translate this to code.

```python
s1 = SchubertCycle(a=[1], k=2, n=4)
s3 = SchubertCycle(a=[3], k=2, n=4)
s1*s1*s1*s3
```

Which yields the output

```
σ_(3,3)
```

Hence, there is a single line that intersects 3 planes and passes through a fixed point in `P^4`.

## Future Work




## Citations

The formulas used in the code come from Sheldon Katz's book on Enumerative Geometry - "Enumerative Geometry and String Theory" which is part of the IAS/Park City Mathematical Subseries of the Student Mathematical Library.
