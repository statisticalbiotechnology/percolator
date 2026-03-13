# Convex q-value spline for PEP estimation

This note describes a q-to-PEP calibration scheme where we fit a smooth,
convex spline to the q-value curve and then recover the local posterior error
probability from its derivative.

## Setup

Let `u` denote the fraction of accepted targets, or equivalently a normalized
rank running from `0` to `1` through the ordered target list. Let

`q(u)`

be the q-value curve over that ordered list, and let

`p(u)`

be the local error rate, i.e. the PEP.

The usual cumulative relationship is

`q(u) = (1 / u) \int_0^u p(t) dt`

for `u > 0`.

Differentiating gives

`u q'(u) + q(u) = p(u)`

so once we have a differentiable estimate of `q(u)`, the implied PEP curve is

`p(u) = q(u) + u q'(u)`.

This is the key identity: q-values are cumulative averages, while PEPs are the
corresponding local rates.

## Why fit a convex q spline

If the local error rate worsens as we move down the ranked list, then `p(u)`
should be nondecreasing in `u`. Under

`p(u) = q(u) + u q'(u)`,

that means the fitted q-curve should be smooth and shaped so the derived local
rate is monotone. A convex spline is a convenient way to enforce that shape:

- `q(u)` stays smooth rather than piecewise constant.
- `q'(u)` is nondecreasing.
- `p(u)` inherits a stable, smooth trend instead of noisy finite differences.

In practice we also clamp both `q(u)` and `p(u)` to `[0, 1]`.

## I-spline parameterization

We want a basis that makes the shape constraints easy to impose. I-splines are
useful because they are monotone basis functions, and nonnegative coefficients
preserve that monotonicity.

One convenient construction is:

1. Represent the slope `q'(u)` with a nonnegative I-spline expansion.
2. Integrate that expansion to obtain `q(u)`.
3. Recover `p(u)` from `q(u) + u q'(u)`.

For example,

`q'(u) = b0 + \sum_j b_j I_j(u),    b0 >= 0, b_j >= 0`

which makes `q'(u)` nondecreasing. Integrating gives

`q(u) = a0 + b0 u + \sum_j b_j \int_0^u I_j(t) dt`.

Then

`p(u) = q(u) + u q'(u)`.

This keeps the optimization linear in the spline coefficients before the final
`p(u)` evaluation, and the nonnegativity constraints are simple box constraints.

## Practical fitting procedure

Given target q-values `q_i` at normalized ranks `u_i`:

1. Sort targets by score from best to worst.
2. Map positions to `u_i` in `(0, 1]`.
3. Build the I-spline design matrix on `u_i`.
4. Fit the convex q-spline coefficients with least squares plus a small ridge
   penalty.
5. Enforce nonnegative spline coefficients to keep the slope monotone.
6. Evaluate both `q(u_i)` and `q'(u_i)`.
7. Compute `p_i = q(u_i) + u_i q'(u_i)`.
8. Clip to `[0, 1]`, and optionally project `p_i` to be nondecreasing if a
   final monotonicity cleanup is wanted.

## Why this is attractive

- It works directly from q-values, so it stays consistent with any upstream
  q-value correction.
- It avoids differencing noisy empirical q-curves.
- The monotonicity/convexity constraints become simple coefficient bounds.
- It matches the existing spline-oriented implementation style in this codebase.

## Implementation notes

- Use normalized `u` rather than raw score values when the goal is specifically
  q-to-PEP conversion.
- Keep an intercept term for the baseline q-level near the top of the list.
- Use a small ridge penalty to stabilize fits when the q-curve has long flat
  stretches.
- The derivative of the spline basis should be evaluated analytically rather
  than by finite differences, since `p(u)` depends directly on `q'(u)`.

The main algorithmic idea is therefore:

fit a convex spline to `q(u)`, then convert it to PEPs through

`p(u) = q(u) + u q'(u)`.


Key identity:
p(u) = q(u) + u q'(u)

### Further implementation constraints:

- q'(u) >= 0
- q''(u) >= 0
- parameterization: q(u)=a + b u + Σ c_j J_j(u)
- b >= 0, c_j >= 0

The codebase already contains:

- ISplineTRRRegressor
- TRR boxed quadratic solver
- Eigen matrices

It makes sense to replace the current I-spline implementation with one implementing this convex parameterization.

