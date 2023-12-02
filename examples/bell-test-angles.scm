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

;;;
;;; Simulation of a two-channel optical Bell test. See, for instance,
;;; https://en.wikipedia.org/w/index.php?title=CHSH_inequality&oldid=1185876217
;;;
;;; The analysis there is wrong, and the experiments run by Aspect et
;;; alia are meaningless. The analysis and simulation here involve no
;;; quantum mechanics and yet results same as those predicted by
;;; quantum mechanics are obtained.
;;;
;;; This is a simple consequence of two facts about mathematics:
;;;
;;;    (1) There is no such thing as a problem that cannot be
;;;        converted to a different subject matter. Thus there is no
;;;        such thing as a problem being ‘quantum’ in character.
;;;
;;;    (2) All methods must reach the same result. Therefore quantum
;;;        mechanics NEVER is necessary to solve a problem.
;;;
;;; All the activity is ‘local realistic’.
;;;
;;; (The alternative to ‘local realistic’ is ‘non-local preposterous’,
;;; which future generations will mock relentlessly. ‘Instantaneous
;;; action at a distance’ had to be and was eliminated both from
;;; Galilean relativity and Newton’s law of universal gravitation, for
;;; otherwise there is a gap in the structure of human knowledge. The
;;; person who did both eliminations was Albert Einstein. But then
;;; ‘instantaneous action at a distance’ was restored by the risible,
;;; incredible incompetence of those who did not take seriously the
;;; argument of Einstein, Podolsky, and Rosen (EPR). A better argument
;;; is (1) and (2) above, but the argument of EPR should have sufficed
;;; until this, likely unshakeable one was found.)
;;;

(import (scheme base)
        (scheme case-lambda)
        (scheme write)
        (only (srfi 144) fl-epsilon)
        (epr-simulation))

(cond-expand
  (chicken (import (format))
           (import (matchable))))

(define bell-test-angles
  `((0 ,π/8)
    (0 ,π3/8)
    (,π/4 ,π/8)
    (,π/4 ,π3/8)))

(define θH 0)
(define θV π/2)

(define (detection-probabilities φ₁ φ₂)
  ;;
  ;; The probabilities of detections are solved for below easily,
  ;; without quantum mechanics. NO OTHER SOLUTIONS ARE POSSIBLE. They
  ;; are required for mathematical consistency.
  ;;
  ;; The calculation via quantum mechanics is merely the following
  ;; made more difficult. One factors the probabilities into complex
  ;; conjugate pairs, then uses those factors as amplitudes of vector
  ;; components in a linear state space. In the end, one multiplies
  ;; conjugate pairs of calculated amplitudes and thus gets back
  ;; probabilities.
  ;;
  ;; The orthodoxy in physics attributes physical meaning to these
  ;; vector components, their linear combinations, and their
  ;; amplitudes. There is no physical meaning. The orthodoxy is
  ;; wrong. They also do not know probability theory nor understand
  ;; what they are doing with their linear algebra.
  ;;
  (define pbs₁ (make-pbs φ₁))
  (define pbs₂ (make-pbs φ₂))
  (let-values (((PH₁V₂ PV₁H₂) (photon-pair-probabilities))
               ((PH₁+ PH₁-) (pbs-probabilities pbs₁ θH))
               ((PV₁+ PV₁-) (pbs-probabilities pbs₁ θV))
               ((PH₂+ PH₂-) (pbs-probabilities pbs₂ θH))
               ((PV₂+ PV₂-) (pbs-probabilities pbs₂ θV)))
    ;; (H + V +) -- horiz (+) at pbs₁  vert (+) at pbs₂
    ;; (H + V -) -- horiz (+) at pbs₁  vert (-) at pbs₂
    ;; etc.
    `(((H + V +) ,(* PH₁V₂ PH₁+ PV₂+))
      ((H + V -) ,(* PH₁V₂ PH₁+ PV₂-))
      ((H - V +) ,(* PH₁V₂ PH₁- PV₂+))
      ((H - V -) ,(* PH₁V₂ PH₁- PV₂-))
      ((V + H +) ,(* PV₁H₂ PV₁+ PH₂+))
      ((V + H -) ,(* PV₁H₂ PV₁+ PH₂-))
      ((V - H +) ,(* PV₁H₂ PV₁- PH₂+))
      ((V - H -) ,(* PV₁H₂ PV₁- PH₂-)))))

(define (photon->symbol phot)
  (if (< (photon-polarization-angle phot) 0.0001) 'H 'V))

(define (simulate-one-event pbs₁ pbs₂)
  (let*-values (((phot₁ phot₂) (photon-pair-source θH θV))
                ((detect₁+ _detect₁-) (pbs-activity pbs₁ phot₁))
                ((detect₂+ _detect₂-) (pbs-activity pbs₂ phot₂)))
    ;; (H + V +) -- horiz (+) at pbs₁  vert (+) at pbs₂
    ;; (H + V -) -- horiz (+) at pbs₁  vert (-) at pbs₂
    ;; etc.
    `(,(photon->symbol phot₁) ,(if detect₁+ '+ '-)
      ,(photon->symbol phot₂) ,(if detect₂+ '+ '-))))

(define *events-per-test-angle* (make-parameter 100000))

(define (detection-frequencies φ₁ φ₂)
  ;; Simulate events and compute frequencies of the different
  ;; detection patterns.
  (define pbs₁ (make-pbs φ₁))
  (define pbs₂ (make-pbs φ₂))
  (define N (*events-per-test-angle*))
  (define NH+V+ 0) (define NH+V- 0)
  (define NH-V+ 0) (define NH-V- 0)
  (define NV+H+ 0) (define NV+H- 0)
  (define NV-H+ 0) (define NV-H- 0)
  (do ((i 0 (+ i 1)))
      ((= i N))
    (match (simulate-one-event pbs₁ pbs₂)
      ('(H + V +) (set! NH+V+ (+ NH+V+ 1)))
      ('(H + V -) (set! NH+V- (+ NH+V- 1)))
      ('(H - V +) (set! NH-V+ (+ NH-V+ 1)))
      ('(H - V -) (set! NH-V- (+ NH-V- 1)))
      ('(V + H +) (set! NV+H+ (+ NV+H+ 1)))
      ('(V + H -) (set! NV+H- (+ NV+H- 1)))
      ('(V - H +) (set! NV-H+ (+ NV-H+ 1)))
      ('(V - H -) (set! NV-H- (+ NV-H- 1)))))
  ;; (H + V +) -- horiz (+) at pbs₁  vert (+) at pbs₂
  ;; (H + V -) -- horiz (+) at pbs₁  vert (-) at pbs₂
  ;; etc.
  `(((H + V +) ,(/ NH+V+ N))
    ((H + V -) ,(/ NH+V- N))
    ((H - V +) ,(/ NH-V+ N))
    ((H - V -) ,(/ NH-V- N))
    ((V + H +) ,(/ NV+H+ N))
    ((V + H -) ,(/ NV+H- N))
    ((V - H +) ,(/ NV-H+ N))
    ((V - H -) ,(/ NV-H- N))))

(define (estimate-correlation detection-freqs)
  ;; Use detection frequencies and trigonometry to estimate the value
  ;; of -cos(2(φ₁-φ₂)=-(cos²(φ₁-φ₂)-sin²(φ₁-φ₂)).
  (define (get-freq pattern)
    (cadr (assoc pattern detection-freqs)))
  (let ((fH+V+ (get-freq '(H + V +)))
        (fH+V- (get-freq '(H + V -)))
        (fH-V+ (get-freq '(H - V +)))
        (fH-V- (get-freq '(H - V -)))
        (fV+H+ (get-freq '(V + H +)))
        (fV+H- (get-freq '(V + H -)))
        (fV-H+ (get-freq '(V - H +)))
        (fV-H- (get-freq '(V - H -))))
    ;; Compute estimates of products of squares of cosines and sines.
    (let ((cos²φ₁sin²φ₂ (+ fH+V+ fV-H-))
          (cos²φ₁cos²φ₂ (+ fH+V- fV-H+))
          (sin²φ₁sin²φ₂ (+ fH-V+ fV+H-))
          (sin²φ₁cos²φ₂ (+ fH-V- fV+H+)))
      ;; Take square roots. All the test angles are in Quadrant I, and
      ;; so only positive square roots will be needed. (Be careful!
      ;; You have to account for the quadrants of φ₁ and φ₂, and so
      ;; sometimes need a NEGATIVE square root when doing this kind of
      ;; calculation.)
      (let ((cosφ₁sinφ₂ (sqrt cos²φ₁sin²φ₂))
            (cosφ₁cosφ₂ (sqrt cos²φ₁cos²φ₂))
            (sinφ₁sinφ₂ (sqrt sin²φ₁sin²φ₂))
            (sinφ₁cosφ₂ (sqrt sin²φ₁cos²φ₂)))
        ;; Use angle-difference identities. See, for instance, the CRC
        ;; Handbook of Mathematical Sciences, 6th edition, page 170.
        (let ((sin<φ₁-φ₂> (- sinφ₁cosφ₂ cosφ₁sinφ₂))
              (cos<φ₁-φ₂> (+ cosφ₁cosφ₂ sinφ₁sinφ₂)))
          ;; That is it. We have everthing we need.
          (- (square sin<φ₁-φ₂>) (square cos<φ₁-φ₂>)))))))

(define (angle->string angle)
  (let* ((angle/π (/ angle π))
         (angle/π*64 (* 64 angle/π))
         (iangle/π*64 (round angle/π*64))
         (diff (abs (- angle/π*64 iangle/π*64)))
         (exact-enough (<= diff (* 500 fl-epsilon
                                   (abs angle/π*64)))))
    (string-append
     "π×" (let ((angle/π (if exact-enough (exact angle/π) angle/π)))
            (number->string angle/π)))))

(define pattern-list
  '((H + V +) (H + V -) (H - V +) (H - V -)
    (V + H +) (V + H -) (V - H +) (V - H -)))

(format #t "~%")

(format #t "~2Tlegend:~%")
(format #t "~2T  (H + V +)  horizontal photon in (+) channel of pbs₁,~%")
(format #t "~2T             vertical photon in (+) channel of pbs₂,~%")
(format #t "~2T  (H + V -)  horizontal photon in (+) channel of pbs₁,~%")
(format #t "~2T             vertical photon in (-) channel of pbs₂, etc.~%")

;;; See
;;; https://en.wikipedia.org/w/index.php?title=CHSH_inequality&oldid=1185876217
;;; Be aware that the contents of the page is mostly wrong. That is
;;; proven by this simulation. This contrast calculation is included
;;; to emphasize our proof.
(define S-nominal 0)   ;; For computing a CHSH contrast.
(define S-estimated 0) ;; For computing a CHSH contrast.
(define i 1)           ;; For computing a CHSH contrast.

(format #t "~%")
(do ((test-angles bell-test-angles (cdr test-angles)))
    ((null? test-angles))
  (let* ((φ₁ (caar test-angles))
         (φ₂ (cadar test-angles))
         (probs-list (detection-probabilities φ₁ φ₂))
         (freqs-list (detection-frequencies φ₁ φ₂))
         (nominal-correlation (- (cos (* 2 (- φ₁ φ₂)))))
         (estimated-correlation (estimate-correlation freqs-list)))
    (set! S-nominal
      ((if (= i 2) - +) S-nominal nominal-correlation))
    (set! S-estimated
      ((if (= i 2) - +) S-estimated estimated-correlation))
    (set! i (+ i 1))
    (format #t "  test angles:  φ₁ = ~A   φ₂ = ~A~%"
            (angle->string φ₁) (angle->string φ₂))
    (format #t "                       nominal     simulated~%")
    (do ((patterns pattern-list (cdr patterns)))
        ((null? patterns))
      (let* ((patt (car patterns))
             (prob (cadr (assoc patt probs-list)))
             (freq (cadr (assoc patt freqs-list))))
        (format #t "  ~A freq~14,5@F~14,5@F~%"
                patt (inexact prob) (inexact freq))))
    (format #t "     correlation~14,5@F~14,5@F~%"
            nominal-correlation estimated-correlation)
    (format #t "~%")))
(format #t "                       nominal     simulated~%")
(format #t "    CHSH S value~14,5@F~14,5@F~%"
        S-nominal S-estimated)
(format #t "~%")
(format #t "  (See https://en.wikipedia.org/w/index.php?title=CHSH_inequality&oldid=1185876217~%")
(format #t "  for a page concerning the claim this result cannot be~%")
(format #t "  gotten without ‘quantum entanglement’. Note ‘quantum~%")
(format #t "  correlation’ is simply the correlation coefficient.~%")
(format #t "  There is nothing ‘quantum’ about it.)~%")
(format #t "~%")
