# epr-simulation

An R⁷RS Scheme toolkit for constructing fully ‘local realistic’
simulations of EPR-B experiments. It is donated to the public domain.

*This toolkit is to help you construct simulations that prove widely
accepted physics wrong, including physics for which the Nobel Prize
was awarded in 2022. Some example simulations are included in the
package.*

## This toolkit versus widespread beliefs

Widespread belief in physics (in 2023) is that experiments of the type
conducted by Alain Aspect are distinctly ‘quantum’ and distinct from
‘classical’ physics. There is supposed to be an ‘entangled pair’
rather than distinct particles, and ‘instantaneous action at a
distance’ resulting in distinct particles.

All of this is easily proven false, for instance empirically by
simulations constructed using this toolkit.

Widespread belief is not evidence in favor of a theory. Belief in
orthodoxy is not scientific method, it is deference to authority. To
conduct experiments by writing simulations is an example of scientific
method.

Another example of scientific method is to make solid logical
arguments, such as the following.

## The impossibility of ‘irreducibility’

That there is no ‘quantum’ physics distinct from physics generally can
actually be shown *logically* by the following argument—

Once a *physics* problem is put into words, it becomes a word
problem in *mathematics*. It is no longer in the realm of
physics. As *mathematics*, the problem is subject to the principles
of *that* field of study. One of those principles is that
mathematics has no fixed subject matter. Any problem in mathematics
can be translated into an endless number of problems concerning
different topics. Furthermore, any method used to solve a problem must
arrive at the same conclusion as any other.

Thus a problem in ‘quantum’ physics can be converted to, say, a
problem about Girl Scouts delivering cookies to their customers, or a
problem about radio signals, or a problem about abstract symbols. The
new problem may be a contrived one, but the conversion *can* be
done. Also, any method used to solve the problem *must* arrive at
the same result as does quantum mechanics, if the solution by quantum
mechanics were itself correct.

The supposed ‘irreducibility’ of ‘quantum’ physics is not merely
false: it is a logically solid impossibility.  That is, there is not
really such a thing as ‘quantum’ physics distinct from ‘classical’
physics: there are, rather, problems that heretofore have mistakenly
been believed soluble ''only'' by quantum mechanical methods, but
which actually can be solved by ''any'' reasonable mathematical
methods. Any method used must reach the same result as any other.

Thus, in the published papers (such as those by J. S. Bell, John
Clauser, etc.) where a ‘local realistic’ solution for EPR-B
correlations is reached that is different from what quantum mechanics
gives, the authors *definitely* have made errors. There is not the
least doubt about it, because the solution from quantum mechanics
*is* correct.

I will prove that below, by deriving the same expression *without*
using quantum mechanics.

## The correlation coefficient of the two-channel optical Bell test

### Photons

I shall model the light as particles, called ‘photons’. Each photon is
plane-polarized at angle 0 or π/2, but we will call the first case H
for horizontal and the latter V for vertical. This notation makes the
mathematics easier to read. The photons are generated in complementary
pairs, polarized thus: (H, V) or (V, H).

Let the left side be subscripted ₁ and the right side subscripted ₂,
in whatever follows.

### Two-channel polarizers

The intensity of the light in a polarizer’s channel, by the Law of
Malus, is inferred to be the probability of a photon being
retransmitted into that channel. Mathematical consistency requires
this inference. The same probability may be deduced from quantum
mechanics.

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
correlation coefficient equal the covariance. That is, no
normalization is needed. It can be computed (with some care) as
follows—

    ρ = ½(+1)(+1) sin²φ₀₁ + ½(+1)(−1) cos²φ₀₁
          + ½(−1)(+1) cos²φ₀₁ + ½(−1)(−1) sin²φ₀₁

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
backwards and found the derivation given above, thus proving that
‘quantum correlation’ is simply the correlation, and that Einstein,
Podolsky, and Rosen were correct.
