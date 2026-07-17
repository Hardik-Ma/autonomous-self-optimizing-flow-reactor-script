# Derivation of the Van de Vusse Physics-Informed Prior Mean Functions

This derives, with every algebraic step shown, the two closed-form expressions used in `piGP_prior_mean.m`:

$$\text{yield}_B = Da\,e^{-Da}, \qquad \text{sel} = \frac{Da\,e^{-Da}}{1-e^{-Da}}$$

---

## 1. Governing equations

Drop the side reaction $2A\to D$ and keep only the consecutive series $A\xrightarrow{k_1}B\xrightarrow{k_2}C$, both first order, isothermal. In a PFR, residence time $\tau$ plays the role of batch time, so the mole balances are:

$$\frac{dc_A}{d\tau} = -k_1 c_A \tag{1}$$
$$\frac{dc_B}{d\tau} = k_1 c_A - k_2 c_B \tag{2}$$
$$\frac{dc_C}{d\tau} = k_2 c_B \tag{3}$$

Initial (inlet) conditions: $c_A(0) = c_{A,0}$, $c_B(0) = 0$, $c_C(0) = 0$.

---

## 2. Solving for $c_A(\tau)$

Equation (1) is separable:
$$\frac{dc_A}{c_A} = -k_1\,d\tau$$

Integrate both sides, left from $c_{A,0}$ to $c_A(\tau)$, right from $0$ to $\tau$:
$$\int_{c_{A,0}}^{c_A(\tau)} \frac{dc_A}{c_A} = -k_1\int_0^\tau d\tau$$
$$\ln c_A(\tau) - \ln c_{A,0} = -k_1\tau$$
$$\ln\!\left(\frac{c_A(\tau)}{c_{A,0}}\right) = -k_1\tau$$

Exponentiate both sides:
$$\boxed{c_A(\tau) = c_{A,0}\,e^{-k_1\tau}} \tag{4}$$

---

## 3. Solving for $c_B(\tau)$ — general case ($k_1 \neq k_2$)

Substitute (4) into (2):
$$\frac{dc_B}{d\tau} = k_1 c_{A,0}e^{-k_1\tau} - k_2 c_B$$

Rearrange into standard linear first-order form $\dfrac{dc_B}{d\tau} + P(\tau)c_B = Q(\tau)$:
$$\frac{dc_B}{d\tau} + k_2 c_B = k_1 c_{A,0}\,e^{-k_1\tau} \tag{5}$$

with $P = k_2$ (constant), $Q(\tau) = k_1 c_{A,0}e^{-k_1\tau}$.

**Integrating factor:**
$$\mu(\tau) = \exp\!\left(\int k_2\,d\tau\right) = e^{k_2\tau}$$

Multiply (5) through by $\mu(\tau)$:
$$e^{k_2\tau}\frac{dc_B}{d\tau} + k_2 e^{k_2\tau}c_B = k_1 c_{A,0}\,e^{k_2\tau}e^{-k_1\tau} \tag{6}$$

**Check the left side is a total derivative** (product rule on $c_B(\tau)e^{k_2\tau}$):
$$\frac{d}{d\tau}\Big[c_B(\tau)e^{k_2\tau}\Big] = \frac{dc_B}{d\tau}e^{k_2\tau} + c_B\cdot k_2 e^{k_2\tau}$$

This is exactly the left side of (6). Good — so (6) becomes:
$$\frac{d}{d\tau}\Big[c_B(\tau)e^{k_2\tau}\Big] = k_1 c_{A,0}\,e^{(k_2-k_1)\tau} \tag{7}$$

(using $e^{k_2\tau}e^{-k_1\tau} = e^{(k_2-k_1)\tau}$)

**Integrate (7) from $0$ to $\tau$** (dummy variable $s$):
$$\int_0^\tau \frac{d}{ds}\Big[c_B(s)e^{k_2 s}\Big]ds = \int_0^\tau k_1 c_{A,0}\,e^{(k_2-k_1)s}\,ds$$

Left side, by the fundamental theorem of calculus:
$$c_B(\tau)e^{k_2\tau} - c_B(0)e^{0} = c_B(\tau)e^{k_2\tau} - 0 = c_B(\tau)e^{k_2\tau}$$

Right side (assuming $k_2\neq k_1$, so this is a genuine exponential integral):
$$k_1 c_{A,0}\left[\frac{e^{(k_2-k_1)s}}{k_2-k_1}\right]_0^\tau = k_1 c_{A,0}\cdot\frac{e^{(k_2-k_1)\tau} - 1}{k_2-k_1}$$

So:
$$c_B(\tau)e^{k_2\tau} = \frac{k_1 c_{A,0}}{k_2-k_1}\Big[e^{(k_2-k_1)\tau} - 1\Big]$$

**Solve for $c_B(\tau)$** by dividing both sides by $e^{k_2\tau}$:
$$c_B(\tau) = \frac{k_1 c_{A,0}}{k_2-k_1}\Big[e^{(k_2-k_1)\tau} - 1\Big]\,e^{-k_2\tau}$$

Distribute $e^{-k_2\tau}$ into the bracket:
$$c_B(\tau) = \frac{k_1 c_{A,0}}{k_2-k_1}\Big[e^{(k_2-k_1)\tau}e^{-k_2\tau} - e^{-k_2\tau}\Big]$$

Simplify the exponent in the first term: $(k_2-k_1)\tau - k_2\tau = k_2\tau - k_1\tau - k_2\tau = -k_1\tau$:
$$\boxed{c_B(\tau) = \frac{k_1 c_{A,0}}{k_2-k_1}\Big[e^{-k_1\tau} - e^{-k_2\tau}\Big]} \qquad (k_1\neq k_2) \tag{8}$$

---

## 4. The degenerate case $k_1 = k_2$

Your locked kinetic set has $k_1 \equiv k_2$ exactly (confirmed numerically earlier — same $k_{10}, E_{a,1}$ as $k_{20}, E_{a,2}$), so (8) is a $0/0$ form. Two independent ways to resolve it — shown both, as a cross-check.

### 4a. Direct method: substitute $k_1=k_2=k$ into the ODE from scratch

Equation (5) becomes, with $k_2\to k$ and using $k_1=k$ in the forcing term too:
$$\frac{dc_B}{d\tau} + k\,c_B = k\,c_{A,0}\,e^{-k\tau}$$

Integrating factor $\mu(\tau) = e^{k\tau}$, multiply through:
$$\frac{d}{d\tau}\Big[c_B(\tau)e^{k\tau}\Big] = k\,c_{A,0}\,e^{k\tau}e^{-k\tau} = k\,c_{A,0}$$

The right side is now a **constant** (the exponentials cancel exactly since the forcing and the decay share the same rate). Integrate from $0$ to $\tau$:
$$c_B(\tau)e^{k\tau} - 0 = k\,c_{A,0}\int_0^\tau ds = k\,c_{A,0}\,\tau$$

Divide by $e^{k\tau}$:
$$\boxed{c_B(\tau) = c_{A,0}\,k\,\tau\,e^{-k\tau}} \tag{9}$$

### 4b. Cross-check: take the limit of the general formula (8) as $k_2\to k_1$

Set $k_2 = k_1+\varepsilon$ and let $\varepsilon\to 0$. The numerator of the bracket in (8):
$$e^{-k_1\tau} - e^{-(k_1+\varepsilon)\tau} = e^{-k_1\tau} - e^{-k_1\tau}e^{-\varepsilon\tau} = e^{-k_1\tau}\Big[1-e^{-\varepsilon\tau}\Big]$$

Expand $e^{-\varepsilon\tau}$ to first order in $\varepsilon$ (small $\varepsilon$): $e^{-\varepsilon\tau} \approx 1-\varepsilon\tau$, so:
$$1-e^{-\varepsilon\tau} \approx \varepsilon\tau$$

Numerator $\approx e^{-k_1\tau}\cdot \varepsilon\tau$. The denominator of (8) is $k_2-k_1=\varepsilon$. So the bracket:
$$\frac{e^{-k_1\tau}-e^{-k_2\tau}}{k_2-k_1} \;\xrightarrow{\varepsilon\to0}\; \frac{e^{-k_1\tau}\cdot\varepsilon\tau}{\varepsilon} = \tau\,e^{-k_1\tau}$$

Substituting back into (8):
$$c_B(\tau) \to k_1 c_{A,0}\cdot \tau\,e^{-k_1\tau}$$

Identical to (9). Both routes agree — good consistency check.

---

## 5. Yield of B

$$\text{yield}_B(\tau) = \frac{c_B(\tau)}{c_{A,0}} = \frac{c_{A,0}\,k\,\tau\,e^{-k\tau}}{c_{A,0}} = k\,\tau\,e^{-k\tau}$$

**Define the Damköhler number** for this first-order step:
$$Da \equiv k\,\tau$$

(dimensionless — note there is *no* concentration term here, unlike a second-order Damköhler number, because the reaction is first order in $c_A$)

Substituting $Da$:
$$\boxed{\text{yield}_B = Da\,e^{-Da}} \tag{10}$$

This is exactly `Da_yield * exp(-Da_yield)` in the code.

---

## 6. Selectivity

### 6a. Mass balance identity (side reaction off)

Add equations (1)+(2)+(3):
$$\frac{d}{d\tau}\big[c_A+c_B+c_C\big] = \underbrace{-k_1c_A}_{(1)} + \underbrace{(k_1c_A - k_2c_B)}_{(2)} + \underbrace{k_2c_B}_{(3)}$$

The right side telescopes to zero:
$$-k_1c_A + k_1c_A - k_2c_B + k_2c_B = 0$$

So $c_A+c_B+c_C$ is constant in $\tau$. Evaluate the constant at $\tau=0$:
$$c_A(0)+c_B(0)+c_C(0) = c_{A,0}+0+0 = c_{A,0}$$

Therefore, for **all** $\tau$:
$$c_A(\tau)+c_B(\tau)+c_C(\tau) = c_{A,0} \quad\Rightarrow\quad c_C(\tau) = c_{A,0}-c_A(\tau)-c_B(\tau) \tag{11}$$

### 6b. Selectivity formula

With the side reaction off, $c_D=0$, so:
$$\text{sel}(\tau) = \frac{c_B(\tau)}{c_B(\tau)+c_C(\tau)+c_D(\tau)} = \frac{c_B(\tau)}{c_B(\tau)+c_C(\tau)}$$

Substitute (11) for $c_C(\tau)$:
$$c_B(\tau)+c_C(\tau) = c_B(\tau) + \big[c_{A,0}-c_A(\tau)-c_B(\tau)\big] = c_{A,0}-c_A(\tau)$$

(the $c_B(\tau)$ terms cancel). So:
$$\text{sel}(\tau) = \frac{c_B(\tau)}{c_{A,0}-c_A(\tau)} \tag{12}$$

### 6c. Substitute the closed forms

Using $c_A(\tau) = c_{A,0}e^{-k\tau}$ from (4):
$$c_{A,0}-c_A(\tau) = c_{A,0}-c_{A,0}e^{-k\tau} = c_{A,0}\big(1-e^{-k\tau}\big)$$

Using $c_B(\tau) = c_{A,0}k\tau e^{-k\tau}$ from (9), substitute both into (12):
$$\text{sel}(\tau) = \frac{c_{A,0}\,k\tau\,e^{-k\tau}}{c_{A,0}\big(1-e^{-k\tau}\big)}$$

Cancel $c_{A,0}$ (appears in both numerator and denominator):
$$\text{sel}(\tau) = \frac{k\tau\,e^{-k\tau}}{1-e^{-k\tau}}$$

Substitute $Da=k\tau$:
$$\boxed{\text{sel} = \frac{Da\,e^{-Da}}{1-e^{-Da}}} \tag{13}$$

This is exactly `series_ratio(Da)` in the code.

---

## 7. Resolving the $Da\to0$ singularity (Taylor expansion, full steps)

Equation (13) is $0/0$ at $Da=0$ (both numerator and denominator vanish). Full series derivation of the removable singularity:

### 7a. Expand $e^{-Da}$

$$e^{-Da} = 1 - Da + \frac{Da^2}{2} - \frac{Da^3}{6} + \frac{Da^4}{24} - \cdots \tag{14}$$

### 7b. Numerator: $N(Da) = Da\,e^{-Da}$

Multiply (14) by $Da$:
$$N(Da) = Da\left(1-Da+\frac{Da^2}{2}-\cdots\right) = Da - Da^2 + \frac{Da^3}{2} - \cdots \tag{15}$$

### 7c. Denominator: $D(Da) = 1-e^{-Da}$

$$D(Da) = 1-\left(1-Da+\frac{Da^2}{2}-\frac{Da^3}{6}+\cdots\right) = Da - \frac{Da^2}{2} + \frac{Da^3}{6} - \cdots$$

Factor out $Da$:
$$D(Da) = Da\left(1-\frac{Da}{2}+\frac{Da^2}{6}-\cdots\right) \tag{16}$$

### 7d. Form the ratio and cancel $Da$

$$f(Da) = \frac{N(Da)}{D(Da)} = \frac{Da\left(1-Da+\frac{Da^2}{2}-\cdots\right)}{Da\left(1-\frac{Da}{2}+\frac{Da^2}{6}-\cdots\right)} = \frac{1-Da+\frac{Da^2}{2}}{1-\frac{Da}{2}+\frac{Da^2}{6}} + O(Da^3)$$

The explicit $Da$ factor cancels top and bottom — this is *why* the singularity is removable (both numerator and denominator vanish at exactly the same rate, $O(Da^1)$).

### 7e. Expand the reciprocal of the denominator

Let $x = \dfrac{Da}{2}-\dfrac{Da^2}{6}$, so the denominator is $1-x$. Using $\dfrac{1}{1-x} \approx 1+x+x^2$ for small $x$:

$$x^2 = \left(\frac{Da}{2}-\frac{Da^2}{6}\right)^2 = \frac{Da^2}{4} - 2\cdot\frac{Da}{2}\cdot\frac{Da^2}{6} + \cdots = \frac{Da^2}{4} + O(Da^3)$$

(the cross term is already $O(Da^3)$, drop it)

$$\frac{1}{1-\frac{Da}{2}+\frac{Da^2}{6}} \approx 1 + \left(\frac{Da}{2}-\frac{Da^2}{6}\right) + \frac{Da^2}{4} + O(Da^3) = 1+\frac{Da}{2}+\left(\frac{1}{4}-\frac{1}{6}\right)Da^2+O(Da^3)$$

Compute the coefficient: $\dfrac{1}{4}-\dfrac{1}{6} = \dfrac{3}{12}-\dfrac{2}{12}=\dfrac{1}{12}$

$$\frac{1}{1-\frac{Da}{2}+\frac{Da^2}{6}} \approx 1+\frac{Da}{2}+\frac{Da^2}{12} + O(Da^3) \tag{17}$$

### 7f. Multiply by the numerator polynomial

$$f(Da) \approx \left(1-Da+\frac{Da^2}{2}\right)\left(1+\frac{Da}{2}+\frac{Da^2}{12}\right)$$

Expand term by term, keeping only powers up to $Da^2$:

| term | product | order |
|---|---|---|
| $1\times1$ | $1$ | $Da^0$ |
| $1\times\frac{Da}{2}$ | $\frac{Da}{2}$ | $Da^1$ |
| $1\times\frac{Da^2}{12}$ | $\frac{Da^2}{12}$ | $Da^2$ |
| $-Da\times1$ | $-Da$ | $Da^1$ |
| $-Da\times\frac{Da}{2}$ | $-\frac{Da^2}{2}$ | $Da^2$ |
| $-Da\times\frac{Da^2}{12}$ | — | $Da^3$, drop |
| $\frac{Da^2}{2}\times1$ | $\frac{Da^2}{2}$ | $Da^2$ |
| $\frac{Da^2}{2}\times\frac{Da}{2}$ | — | $Da^3$, drop |

**Collect by order:**

$Da^0$: $\quad 1$

$Da^1$: $\quad \dfrac{1}{2}-1 = -\dfrac{1}{2}$

$Da^2$: $\quad \dfrac{1}{12}-\dfrac{1}{2}+\dfrac{1}{2} = \dfrac{1}{12}$ (the $-\frac12$ and $+\frac12$ cancel exactly, leaving just $\frac{1}{12}$)

$$\boxed{f(Da) \approx 1 - \frac{Da}{2} + \frac{Da^2}{12} + O(Da^3)} \tag{18}$$

This confirms $\text{sel}\to1$ as $Da\to0$ — physically correct, since at vanishing residence time whatever trace of $A$ has reacted is still essentially all $B$ (no time for the second step to consume it). This is exactly the branch coded for $Da<10^{-4}$ in `series_ratio(Da)`.

---

## 8. Reintroducing temperature (Arrhenius)

$k$ in $Da=k\tau$ is temperature-dependent. Reference-centred Arrhenius form (same convention as `f_van_de_vusse_kinetics.m`):
$$k(T) = k_{\text{ref}}\exp\!\left[-\frac{E_a}{R}\left(\frac{1}{T}-\frac{1}{T_{\text{ref}}}\right)\right]$$

so:
$$Da(\tau,T) = k_{\text{ref}}\exp\!\left[-\frac{E_a}{R}\left(\frac1T-\frac1{T_{\text{ref}}}\right)\right]\tau_{\text{sec}} \tag{19}$$

Substituting (19) into (10) and (13) gives the final forms used in the prior, before the learned $\alpha,\beta$ correction:
$$\text{yield}_{\text{phys}}(\tau,T) = \alpha_y\big[Da\,e^{-Da}\big]+\beta_y, \qquad \text{sel}_{\text{phys}}(\tau,T) = \alpha_s\left[\frac{Da\,e^{-Da}}{1-e^{-Da}}\right]+\beta_s$$

with $Da$ evaluated separately for each objective, using its own $k_{\text{ref}},E_a$ from $\theta_{\text{yield}}$ or $\theta_{\text{sel}}$.

---

## Summary of what's exact vs. approximate

| Step | Status |
|---|---|
| $c_A(\tau)=c_{A,0}e^{-k\tau}$ | Exact, given side reaction off |
| $c_B(\tau)=c_{A,0}k\tau e^{-k\tau}$ | Exact, given side reaction off *and* $k_1=k_2$ (both true for your locked kinetics) |
| $\text{yield}_B=Da\,e^{-Da}$ | Exact algebraic consequence of the above |
| $\text{sel}=Da\,e^{-Da}/(1-e^{-Da})$ | Exact algebraic consequence of the above |
| Taylor branch near $Da=0$ | Exact local approximation of an exact formula (not approximating the physics, just avoiding floating-point $0/0$) |
| $\alpha,\beta$ correction | The only place approximation enters — absorbing the omitted $2A\to D$ side reaction |
