#| -*- encoding: utf-8; -*-

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
|#

(import (scheme base)
        (scheme write)
        (epr-simulation))

(cond-expand
  (chicken (import (format))
           (import (matchable))))

(define bell-test-angles
  `((0 π/8)
    (0 π3/8)
    (π/4 π/8)
    (π/4 π3/8)))

(define θ₁ 0)
(define θ₂ π/2)

(define √2 (sqrt 2))

(define (quantum-mechanical-probabilities φ₁ φ₂)
  ;;
  ;; These calculations are just ordinary probability calculations,
  ;; but conducted with square roots of probabilities and then squared
  ;; in a final step. In quantum mechanics one would write these steps
  ;; as operations on vectors in a linear state space, but that makes
  ;; no difference. Physics orthodoxy has bestowed physical meaning on
  ;; the vectors in this state space, but we see here that the vectors
  ;; are simply a bookkeeping method for probability calculations.
  ;;
  ;; The ‘physical meaning’ of a superposition is not worth a moment’s
  ;; contemplation, and certainly is not ‘entanglement’. The
  ;; superposition is merely a way of grouping together some numbers
  ;; for calculational convenience.
  ;;
  ;; We could have done very similar, more straightforward
  ;; calculations using photon-pair-probabilities and
  ;; pbs-probabilities.
  ;;
  (define pbs₁ (make-pbs φ₁))
  (define pbs₂ (make-pbs φ₂))
  (let*-values (((Aphot₁ Aphot₂) (photon-pair-amplitudes))
                ((Apbs₁₁+ Apbs₁₁-) (pbs-amplitudes pbs₁ θ₁))
                ((Apbs₁₂+ Apbs₁₂-) (pbs-amplitudes pbs₁ θ₂))
                ((Apbs₂₁+ Apbs₂₁-) (pbs-amplitudes pbs₂ θ₁))
                ((Apbs₂₂+ Apbs₂₂-) (pbs-amplitudes pbs₂ θ₂)))
    `((,φ₁ ,θ₁ + ,(square (* √2 Aphot₁ Apbs₁₁+)))
      (,φ₁ ,θ₁ - ,(square (* √2 Aphot₁ Apbs₁₁-)))
      (,φ₁ ,θ₂ + ,(square (* √2 Aphot₂ Apbs₁₂+)))
      (,φ₁ ,θ₂ - ,(square (* √2 Aphot₂ Apbs₁₂-)))
      (,φ₂ ,θ₁ + ,(square (* √2 Aphot₁ Apbs₂₁+)))
      (,φ₂ ,θ₁ - ,(square (* √2 Aphot₁ Apbs₂₁-)))
      (,φ₂ ,θ₂ + ,(square (* √2 Aphot₂ Apbs₂₂+)))
      (,φ₂ ,θ₂ - ,(square (* √2 Aphot₂ Apbs₂₂-))))))

(define (quantum-mechanical-correlation probabilities)
  ;;
  ;; The correlation coefficient, assigning +1 to (+) detections and
  ;; -1 to (-) detections. With this assignment, the correlation
  ;; coefficient equals the covariance.
  ;;
  ;; By the way, one could also use the closed solution
  ;; -cos(2(θ₁-θ₂))=-(cos²(θ₁-θ₂)-sin²(θ₁-θ₂)).
  ;;
  (match probabilities
    (((_ _ _ P11+) (_ _ _ P11-) (_ _ _ P12+) (_ _ _ P12-)
      (_ _ _ P21+) (_ _ _ P21-) (_ _ _ P22+) (_ _ _ P22-))
     (let ((s1122++ (* 1/2 P11+ P22+))
           (s1122+- (* 1/2 P11+ P22-))
           (s1122-+ (* 1/2 P11- P22+))
           (s1122-- (* 1/2 P11- P22-))
           (s1221++ (* 1/2 P12+ P21+))
           (s1221+- (* 1/2 P12+ P21-))
           (s1221-+ (* 1/2 P12- P21+))
           (s1121-- (* 1/2 P12- P21-)))
       (let ((sum-alike (+ s1122++ s1122-- s1221++ s1121--))
             (sum-unalike (+ s1122+- s1122-+ s1221+- s1221-+)))
         (- sum-alike sum-unalike))))))

(define (simulate-one-event φ₁ φ₂)
  (define pbs₁ (make-pbs φ₁))
  (define pbs₂ (make-pbs φ₂))
  (let*-values (((phot₁ phot₂) (photon-pair-source θ₁ θ₂))
                ((detect₁+ _detect₁-) (pbs-activity pbs₁ phot₁))
                ((detect₂+ _detect₂-) (pbs-activity pbs₂ phot₂)))
    `((,phot₁ ,(if detect₁+ '+ '-))
      (,phot₂ ,(if detect₂+ '+ '-)))))

(write (quantum-mechanical-probabilities 0 π/8))(newline)
(write (quantum-mechanical-correlation (quantum-mechanical-probabilities 0 π/8)))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
