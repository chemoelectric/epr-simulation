# epr-simulation

An R⁷RS Scheme toolkit for constructing fully ‘local realistic’
simulations of EPR-B experiments. It is donated to the public domain.

*This toolkit is to help you construct simulations that prove widely
accepted physics wrong, including physics for which the Nobel Prize
was awarded in 2022. Some example simulations are included in the
package.*

## Scheme implementations currently supported

- [CHICKEN 5](http://call-cc.org/)

## This toolkit versus widespread beliefs

Widespread belief in physics (in 2023) is that [experiments of the
type conducted by Alain
Aspect](https://en.wikipedia.org/w/index.php?title=Aspect%27s_experiment&oldid=1187336241)
are distinctly ‘quantum’ and distinct from ‘classical’ physics. There
is supposed to be an ‘entangled pair’ rather than distinct particles,
and ‘instantaneous action at a distance’ resulting in distinct
particles.

All of this is easily proven false, for instance empirically by
simulations constructed using this toolkit.

Widespread belief is not evidence in favor of a theory. Belief in
orthodoxy is not scientific method, it is deference to authority. To
conduct experiments by writing simulations is an example of scientific
method.

Another example of scientific method is to make solid logical
arguments, such as the following.

## The impossibility of ‘quantum’ physics being distinct

That there is no ‘quantum’ physics distinct from physics generally can
actually be shown *logically* by the following argument—

Once a *physics* problem is put into words, it becomes a word problem
in *mathematics*. It is no longer in the realm of physics. As
*mathematics*, the problem is subject to the principles of *that*
field of study. One of those principles is that mathematics has no
fixed subject matter. Any word problem in mathematics can be
translated into an endless number of problems concerning different
topics. Furthermore, any method used to solve a problem must arrive at
the same conclusion as any other.

Thus a problem in ‘quantum’ physics can be converted to, say, a
problem about Girl Scouts delivering cookies to their customers, or a
problem about radio signals, or a problem about abstract symbols. The
new problem may be a contrived one, but the conversion *can* be done:
for example, a person can devise a procedure to replace any ‘photon’
in a word problem with a ‘frog’ and any ‘polarizer’ with a ‘frog
channeler’. Even such a trivial change removes the ‘quantum’ from a
‘quantum’ physics problem. Also, any method used to solve the problem
*must* arrive at the same result as does quantum mechanics, if the
solution by quantum mechanics were itself correct. Otherwise
mathematics would be uselessly inconsistent, and we know from
experience that it is not.

What is more, quantum mechanics is itself not a fundamental
mathematical field of study, but an application of linear algebra,
probability theory, etc. Any problem in ‘quantum’ physics surely could
be solved by someone aware of those fields of study but not of quantum
mechanics itself.

Therefore the supposed distinction of ‘quantum’ physics problems from
‘classical’ physics problems is not merely false: it is so implausible
we might as well declare it impossible. *There does not exist
‘quantum’ physics distinct from ‘classical’ physics*. The distinction
is a myth. Rather, there are problems that have been believed soluble
*only* by quantum mechanical methods, but which actually can be solved
by *any* reasonable mathematical methods—*all of which must reach the
same result*.

Thus, in the published papers (such as those by J. S. Bell, John
Clauser, etc.) where a ‘local realistic’ solution for EPR-B
correlations is reached that is different from what quantum mechanics
gives, the authors *definitely* have made errors, because the solution
from quantum mechanics *is* correct.

I will prove that below, by deriving the same expression
independently.

## The correlation coefficient of the two-channel optical Bell test

(*Note: I use the name* ρ *for the correlation coefficient because
that is a common name for it in probability theory. This usage has no
connection to the use of* ρ *in physics to represent a density
function.*)

### [Photons](https://en.wikipedia.org/w/index.php?title=Photon&oldid=1188018882)

I shall model the light as particles, called ‘photons’. Each photon is
plane-polarized at angle 0 or π/2, but we will call the first case H
for horizontal and the latter V for vertical. This notation makes the
mathematics easier to read. The photons are generated in complementary
pairs, polarized thus: (H, V) or (V, H).

Let the left side be subscripted ₁ and the right side subscripted ₂,
in whatever follows.

### [Two-channel polarizers](https://en.wikipedia.org/w/index.php?title=Polarizer&oldid=1173590658#Beam-splitting_polarizers)

The intensity of the light in a polarizer’s channel, by the [Law of
Malus](https://en.wikipedia.org/w/index.php?title=Polarizer&oldid=1173590658#Malus's_law_and_other_properties),
is inferred to be the probability of a photon being retransmitted into
that channel. Mathematical consistency requires this inference. The
same probability may be deduced from quantum mechanics.

Let the angle setting of the polarizer by φ. Then the probability of a
photon being retransmitted into the (+) channel is cos²(φ) for an H
photon or sin²(φ) for a V photon, and into the (−) channel is sin²(φ)
for an H photon or cos²(φ) for a V photon.

### Removing one of the polarizers from the problem

Let us set the left and right side polarizers respectively to φ₀₁ ∈ ℜ
and φ₀₂ = 0. In that case on the right side every H photon is
retransmitted into the (+) channel and every V photon into the
(−) channel. Thus we can ignore that there is a polarizing beam
splitter on that side, and simply map

    H ↦ (+)
    V ↦ (−)

Assign the value +1 to a photon detected in a (+) channel and −1 to a
photon detected in a (−) channel. This choice of values makes the
[correlation coefficient equal the
covariance](https://en.wikipedia.org/w/index.php?title=Covariance_and_correlation&oldid=1144835290). That
is, no normalization is needed. It can be computed (with some care) as
follows—

    ρ = ½(+1)(+1)sin²φ₀₁ + ½(+1)(−1)cos²φ₀₁
          + ½(−1)(+1)cos²φ₀₁ + ½(−1)(−1)sin²φ₀₁

Using a double-angle identity, and taking advantage of φ₀₂ = 0, this
becomes

    ρ = −(cos²φ₀₁ − sin²φ₀₁)
      = −cos(2φ₀₁)
      = −cos(2(φ₀₁ − φ₀₂))

In the last expression, φ₀₁ − φ₀₂ is the difference between the two
polarizer settings. It would appear that the correlation coefficient
is invariant as long as this difference is maintained. *This
invariance is profound. It is what gives us the solution without
quantum mechanics.*

Let us make the invariance more formal. Let Δφ ∈ ℜ, φ₁ = φ₀₁ + Δφ, and
φ₂ = φ₀₂ + Δφ = Δφ.  Clearly φ₁ and φ₂ can be *any* settings of the
polarizers: φ₁ − φ₂ = (φ₀₁ +  Δφ) − (φ₀₂ +  Δφ) = φ₀₁ − φ₀₂ is their
difference, and Δφ is a shared rotational offset.

### The general solution, derived without quantum mechanics

Thus the general formula for the correlation coefficient is

    ρ = −cos(2(φ₁ − φ₂))

for any φ₁, φ₂ ∈ ℜ.

This result is the same as that obtained by quantum mechanics. It had
to be: all mathematical methods, applied to the same problem, must
reach the same result.

Physicists currently (2023) compute this function *only* by quantum
mechanics, call it ‘quantum correlation’, and think it peculiarly
‘quantum’. However, let us suppose they had noticed you could rotate
the polarizer settings together, without affecting the correlation, so
one of the polarizers was at zero degrees. One could then have worked
backwards and found the derivation given above—thus proving ‘quantum
correlation’ is simply the correlation, and that [Einstein, Podolsky,
and Rosen](https://doi.org/10.1103/PhysRev.47.777) were correct:
‘quantum’ physics is *ordinary physics* that has been obfuscated.

Physicists finding the derivation above would have discovered, as
well, that the fabled [‘hidden
variables’](https://en.wikipedia.org/w/index.php?title=Hidden-variable_theory&oldid=1188516029)
are a red herring. One should erase the phrase from the lexicon. It
is fallacious to assume that, because quantum mechanics deals with
probabilities and statistical mechanics deals with probabilities, then
quantum mechanics must reduce to statistical mechanics with ‘hidden
variables’. In Aspect-style experiments, quantum mechanics reduces to
simple division of outcomes according to the Law of Malus.

## Toolkit reference

### Constant values

    [constant] π/180 π/8 π/4 π3/8 π/2 π3/4 π

Approximations of fractions of π, stored as the values of variables.

### Angular conversions

    [procedure] (radians->degrees angle-in-radians)
    [procedure] (degrees->radians angle-in-degrees)

The toolkit generally measures angles in radians. These procedures can
help you convert to and from degrees.

    [procedure] (radians->string angle-in-radians)
    [procedure] (radians->string angle-in-radians tolerance)

Converts an angle in radians to a string such as `"π×3/8"` or
`"π×-0.1234"`. If `tolerance` is specified, it represents a multiple
of one half of `fl-epsilon`, and determines whether a number is
rounded off to an exact fraction. The value of `tolerance` defaults to
`1000`.

    [procedure] (string->radians string)

Converts a string such as `"π×3/8"`, `"π3/8"`, `"π×-0.1234"`, or just
`"0.1234"` to a number representing an angle in radians.

### Numerical utilities

    [procedure] (approx= x y)
    [procedure] (approx= x y tolerance)

Compare two real numbers `x` and `y` for approximate equality. If
`tolerance` is specified, it represents a multiple of one half of
`fl-epsilon`. The value of `tolerance` defaults to `1000`.

### Tensors

In the context of this toolkit, *tensors* are abstract vectors that
are ordered tuples of abstract vectors, and the basis tensors are
ordered tuples of basis vectors. The ordered tuples of basis vectors
will be represented as strings, with commas as separators of the parts
of the tuples. A ‘length’ is represented as the `car` of a pair and
the corresponding basis tensor is the `cdr` of the pair. A linear
combination of such primitive tensors is represented by simply putting
them in a list. The order of the list is immaterial.

Here, for example, is a tensor that, in the notations of quantum
mechanics, represents a bit that with equal probability is either `0`
or `1`:

    `((,(sqrt 0.5) . "0") (,(sqrt 0.5) . "1"))

Here is a three-bit register that is either `100` or `011`, with
probabilities 0.75 and 0.25, respectively:

    `((,(sqrt 0.75) . "1,0,0") (,(sqrt 0.25) . "0,1,1"))

The tuples could be ordered differently, if one wished (and of course
the list order is immaterial):

    `((,(sqrt 0.25) . "1,1,0") (,(sqrt 0.75) . "0,0,1"))

In general the ‘lengths’ of tensor components are allowed to be
complex numbers.

Often, in quantum mechanics, the ‘lengths’ are called
*amplitudes*. When multiplied by its complex conjugate, such an
*amplitude* yields a probability. Storing complex conjugate roots of
probabilities as ‘lengths’ of vector components becomes a form of
bookkeeping. Quantum mechanics can end up as nothing more than an
obfuscation of probability theory.

This is, in fact, what has happened with EPR-B experiments. Thus the
toolkit provides means for comparing probability theory calculations
with quantum mechanical calculations, and showing that they arrive at
the same results.

