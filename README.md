# Pieri's

This is a tool for computing the cohomology ring structure of the Grassmannian G(2,n+1) using Schubert calculus. Namely, the tool can compute the cohomology groups using Schubert classes which serve as generating cocycles in the cohomology ring. Using Pieri's formula, we can compute the product between Schubert classes.
Currently, due to the implementation of Schubert cycles and Pieri's formula, multiplication only works when one of the multiplicands is a class of the form `SchubertCalculus(a=a,b=0,n=n)`; although, this will be fixed shortly. 

## Documentation

For some basic examples of usage:

```python
# To define a Schubert cycle, we use the SchubertCycle class
s1 = SchubertCycle(a=1,b=0,n=3)
# We can then multiply cycles together or add them using * and +
s1*s1
σ_(1,1)+σ_2 # output of product
s1 + s1
2σ_1 # output of sum
```

We can also find the generators for the cohomology ring

```python
compute_cohomology_groups(n=3)
[[1], [], [σ_1], [], [σ_2, σ_(1,1)], [], [σ_(2,1)], [], [σ_(2,2)]]
# The output of this function is stored as a list of lists where the i-th list are the generators of the i-th cohomology class
```

While this can be read easily from the output above, if we want a prettified version of the cohomology groups with generators, we can do this easily

```python
groups = compute_cohomology_groups(n=3)
print_cohomology_groups(groups)
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

Finally, we can compute the multiplication table (although the accuracy of the computation is poor for n > 3; namely, when the Pieri's formula implementation does not work, and we can not take advantage of Poincaré duality)

```python
groups = compute_cohomology_groups(n=3)
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

## Applications

Finally, let's use the software to solve a couple enumerative geometry problems.

### Problem 1

Q: How many lines pass through the 4 distinct lines in projective 3-space, P^3?

To solve this, we first translate the problem to Schubert calculus. The set of lines that pass through a given line is the Schubert cycle σ_1. We can then take the intersection pairing on σ_1^4 to find the number of lines passing through any four given lines (in general position).
Let's translate this to code.

```python
s1 = SchubertCycle(a=1,b=0,n=3)
s1**4
# Note that this exponential works because we can use Poincaré duality to compute (s1*s1)*(s1*s1), in general, the following line is equivalent but works in more situations
((s1*s1)*s1)*s1
```

In either case, we obtain the output

```
2σ_(2,2)
```

And the intersection pairing simply yields the scalar 2. Hence, there are 2 lines that intersect 4 arbitrary lines in P^3

### Problem 2
Q: How many lines pass through 6 distinct planes in projective 4-space, P^4?

Equivalently the problem boils down to computing σ_1^6 in Schubert calculus since σ_1 can be though of as the set of lines passing through a particular plane.
Let's translate this to code.

```python
s1 = SchubertCycle(a=1,b=0,n=4)
((((s1*s1)*s1)*s1)*s1)*s1
```

Which yields the output

```
5σ_(3,3)
```

The intersection pairing gives the scalar 5. Hence, there are 5 lines that intersect 6 arbitrary planes in P^4.

### Problem 3
Q: How many lines pass through 3 distinct planes and a specific point in projective 4-space, P^4?

Similar to problem 2, the problem boils down to computing (σ_1^3)*(σ_3). We can translate this to code.

```python
s1 = SchubertCycle(a=1,b=0,n=4)
s3 = SchubertCycle(a=3,b=0,n=4)
((s1*s1)*s1)*s3
```

Which yields the output

```
σ_(3,3)
```

Hence, there is a single line that intersects 3 planes and passes through a fixed point in P^4.

## Future Work

I would like to generalize this implementation for a general Grassmannian G(k,n+1), and of course, fix the multiplication errors listed above.


## Citations

The formulas used in the code come from Sheldon Katz's book on Enumerative Geometry - "Enumerative Geometry and String Theory" which is part of the IAS/Park City Mathematical Subseries of the Student Mathematical Library.
