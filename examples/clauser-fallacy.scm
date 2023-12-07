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
;;;               Author : Barry Schwartz
;;; Date first completed : 3 December 2023
;;;
;;; Simulation of what the mistaken analysis of John Clauser and
;;; others might ACTUALLY represent. See, for instance,
;;; https://en.wikipedia.org/w/index.php?title=CHSH_inequality&oldid=1185876217
;;;
;;; Let us examine J. S. Bell’s 1971 derivation.
;;;
;;; First, let us dispose of λ. There is no ‘hidden variable’ and
;;; there is no need for statistical mechanics or a state space. Even
;;; if Albert Einstein and friends thought there might be such
;;; complicated mechanisms underneath, there are not. ‘Hidden
;;; variables’ turns out to be a red herring, and state space an
;;; overkill. One can use ordinary probability theory, without ever
;;; having studied statistical mechanics.
;;;
;;; (Bell uses the symbol λ similarly in his famous ‘Bertlmann’s
;;; Socks’ lecture, to distract us from his failure to provide a lemma
;;; proving that P(x)=P(x|y)≜P(x∧y)/P(y) for certain propositions x
;;; and y.)
;;;
;;; More importantly, Bell’s expression for E(a,b) (or P(a,b) in
;;; Bell’s original notation) surely is NOT an expectation for an
;;; EPR-B experiment.
;;;
;;;    IT IS NOT A CORRECT EXPECTATION UNLESS THE PROBABILITY DENSITY
;;;    FUNCTION IS CORRECT.
;;;
;;; Whatever Bell is using surely is not something he devised from
;;; knowledge of probability theory, but something inappropriate he
;;; grabbed from his knowledge of statistical mechanics, etc. PROBABLY
;;; HE IS ACTUALLY MODELING A SINGLE PHOTON PASSING THROUGH TWO
;;; POLARIZERS, or some similar situation.
;;;
;;; This should be intuitive once pointed out. It would explain why,
;;; in the analyses of Clauser, the correlations are supposed to
;;; disappear when a Bell test angle is changed from zero to π/4. An
;;; ideal polarizing beam splitter set to zero could literally be left
;;; out of the apparatus, without effect. One will get the same
;;; results as in the EPR-B experiment. A polarizer set to π/4, on the
;;; other hand, will obstruct the normal flow of the particle
;;; beam. Thus, according to the INCORRECT analysis, the correlations
;;; vanish.
;;;
;;; Before the author of this program found how to derive the
;;; correlation coefficient of the Bell test CORRECTLY, by using a
;;; trivial transformation of the argument angles, he found how to do
;;; it non-trivially, using a joint probability density function
;;; (pdf). His method was complicated and EXTREMELY error prone. There
;;; were always little bugs to stamp out, but probably they could all
;;; have been stamped out. The pdf represented the rotational
;;; invariance as a delta function and probably looked nothing like
;;; what Bell had in mind. My point is that one cannot simply multiply
;;; a bunch of +1 and -1, using any probability density one finds in a
;;; physics textbook, written in a notation not used by anyone but
;;; physicists, and expect to get correct results. It is necessary to
;;; know what one is doing, and that requires knowledge of probability
;;; theory.
;;;
;;; The derivation using a simple transformation of the angles is so
;;; much better that one ought never again talk about doing it my old
;;; way or Bell’s way. I myself try to avoid too much exposure to
;;; incorrect analyses of the EPR-B experiment.
;;;
;;;      *  *  *
;;;
;;; The following simulation is of one possible arrangement of two
;;; polarizers in series, and it gives a CHSH ‘contrast’ |S|=√2 < 2.
;;;
;;; Also included is an example of faulty correlation calculation. The
;;; author found himself writing it after some exposure to an
;;; incorrect analysis of the EPR-B experiment. In horror he realized
;;; what he had done. However, he saved the mistake for inclusion in
;;; this example program.
;;;

(import (scheme base)
        (scheme write)
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
  ;; The probabilities of detections are computed below without
  ;; quantum mechanics. This is the solution demanded by probability
  ;; theory. It is merely the division of possible events into their
  ;; proper proportions.
  ;;
  ;; The calculation via quantum mechanics is merely the following,
  ;; obfuscated. The orthodoxy in physics attributes physical meaning
  ;; to this obfuscated, more difficult calculation. There is no
  ;; physical meaning. It is simply a peculiar way to do the same
  ;; calculation as here.
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
    (let ((probs `(((H + V +) ,(* PH₁V₂ PH₁+ PV₂+))
                   ((H + V -) ,(* PH₁V₂ PH₁+ PV₂-))
                   ((H - V +) ,(* PH₁V₂ PH₁- PV₂+))
                   ((H - V -) ,(* PH₁V₂ PH₁- PV₂-))
                   ((V + H +) ,(* PV₁H₂ PV₁+ PH₂+))
                   ((V + H -) ,(* PV₁H₂ PV₁+ PH₂-))
                   ((V - H +) ,(* PV₁H₂ PV₁- PH₂+))
                   ((V - H -) ,(* PV₁H₂ PV₁- PH₂-)))))

      ;; Sanity check: verify that the probabilities add up to one.
      (check-probabilities (map cadr probs))

      probs)))

(define (compute-correlation probabilities)
  ;;
  ;; The correlation coefficient, assigning +1 to (+) detections and
  ;; -1 to (-) detections. With this assignment, the correlation
  ;; coefficient equals the covariance.
  ;;
  ;; This will be a ‘Clauser-style’ calculation, for a single photon
  ;; passing through two polarizing beam splitters.
  ;;
  ;; (To get the correct answer for the Bell test, one CAN use the
  ;; following calculation, but RESTRICTED TO φ₂=0 and using φ₁ to
  ;; represent the difference between settings. This is then the
  ;; general answer, believe it or not. Adding any necessary number Δφ
  ;; to both φ₁ and φ₂ covers all possible cases. This should have
  ;; become obvious from the correlation expression itself
  ;; [-cos(2(φ₁-φ₂))] and worked backwards to a proof of EPR! The
  ;; expression is invariant under such rotations by Δφ.)
  ;;
  (define (get-prob pattern)
    (cadr (assoc pattern probabilities)))
  (let ((PH+V+ (get-prob '(H + V +)))
        (PH+V- (get-prob '(H + V -)))
        (PH-V+ (get-prob '(H - V +)))
        (PH-V- (get-prob '(H - V -)))
        (PV+H+ (get-prob '(V + H +)))
        (PV+H- (get-prob '(V + H -)))
        (PV-H+ (get-prob '(V - H +)))
        (PV-H- (get-prob '(V - H -))))
    (- (+ PH+V+ PH-V- PV+H+ PV-H-)
       (+ PH+V- PH-V+ PV+H- PV-H+))))

(define (photon->symbol phot)
  (if (< (photon-polarization-angle phot) 0.0001) 'H 'V))

(define (simulate-one-event pbs₁ pbs₂)
  (let*-values (((phot₁ phot₂) (photon-pair-source θH θV))
                ;; phot₁ simply gets counted. (The following is
                ;; equivalent to a PBS set to zero.)
                ((detect₁+) (eq? (photon->symbol phot₁) 'H))
                ;; But phot₂ goes through two PBSes. We will simply
                ;; pass the photon from either channel of pbs₁ to the
                ;; input of pbs₂, without change of orientation.
                ((phot₂+ phot₂-) (pbs-activity pbs₁ phot₂))
                ((detect₂+ _detect₂-)
                 (if phot₂+
                     (pbs-activity pbs₂ phot₂+)
                     (pbs-activity pbs₂ phot₂-))))
    ;; (H + V +) -- horiz (+) at pbs₁  vert (+) at pbs₂
    ;; (H + V -) -- horiz (+) at pbs₁  vert (-) at pbs₂
    ;; etc.
    `(,(photon->symbol phot₁) ,(if detect₁+ '+ '-)
      ,(photon->symbol phot₂) ,(if detect₂+ '+ '-))))

(define *events-per-test-angle* (make-parameter 1000000))

(define (detection-frequencies φ₁ φ₂)
  ;; Simulate events and compute frequencies of the different
  ;; detection patterns.

  ;; pbs₁ outputs a photon that we must represent as a <photon>
  ;; record. In this simulation, the photon retransmitted by pbs₁ will
  ;; be polarized according to the setting of pbs₁, and we will assume
  ;; pbs₂ is lined up with pbs₁.
  (define pbs₁ (make-pbs φ₁ #t #t))

;;;;; FIXME: Some other things to try.
;;;;;        This, for instance, also gives |S|=sqrt(2).
  ;; (define change-pola
  ;;    (make-photon-polarization-angle-changer
  ;;     (lambda (θ) (- π/2 θ))))
  ;; (define pbs₁ (make-pbs φ₁ (lambda (_ phot) (change-pola phot))
  ;;                        (lambda (_ phot) (change-pola phot))))

  ;; pbs₂ outputs photons into photodetectors and so can return
  ;; booleans instead of <photon> records.
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

(define (estimate-correlation φ₁ φ₂ detection-freqs)
  ;; Use detection frequencies to estimate the value of -cos(2(φ₁-φ₂).
  (define (get-freq pattern)
    (cadr (assoc pattern detection-freqs)))
  (estimate-pair-correlation
   φ₁ φ₂ 'optical 'complementary
   `(,(get-freq '(H + V +)) ,(get-freq '(H + V -))
     ,(get-freq '(H - V +)) ,(get-freq '(H - V -))
     ,(get-freq '(V + H +)) ,(get-freq '(V + H -))
     ,(get-freq '(V - H +)) ,(get-freq '(V - H -)))))

(define pattern-list
  '((H + V +) (H + V -) (H - V +) (H - V -)
    (V + H +) (V + H -) (V - H +) (V - H -)))

(format #t "~%")
(format #t "  legend:~%")
(format #t "    (H + V +)  horizontal photon in (+) channel of pbs₁,~%")
(format #t "               vertical photon in (+) channel of pbs₂,~%")
(format #t "    (H + V -)  horizontal photon in (+) channel of pbs₁,~%")
(format #t "               vertical photon in (-) channel of pbs₂, etc.~%")

;;; See
;;; https://en.wikipedia.org/w/index.php?title=CHSH_inequality&oldid=1185876217
(define S-estimated 0) ;; For computing a CHSH contrast.
(define i 1)           ;; For computing a CHSH contrast.

(format #t "~%")
(do ((test-angles bell-test-angles (cdr test-angles)))
    ((null? test-angles))
  (let* ((φ₁ (caar test-angles))
         (φ₂ (cadar test-angles))
         (freqs-list (detection-frequencies φ₁ φ₂))
         (estimated-correlation
          (estimate-correlation φ₁ φ₂ freqs-list)))
    (set! S-estimated
      ((if (= i 2) - +) S-estimated estimated-correlation))
    (set! i (+ i 1))
    (format #t "  test angles:  φ₁ = ~A   φ₂ = ~A~%"
            (radians->string φ₁) (radians->string φ₂))
    (do ((patterns pattern-list (cdr patterns)))
        ((null? patterns))
      (let* ((patt (car patterns))
             (freq (cadr (assoc patt freqs-list))))
        (format #t "  ~A freq~14,5@F~%" patt (inexact freq))))
    (format #t "     correlation~14,5@F~%" estimated-correlation)
    (format #t "~%")))
(format #t "    CHSH S value~14,5@F~%" S-estimated)
(format #t "~%")
(format #t "  See https://en.wikipedia.org/w/index.php?title=CHSH_inequality&oldid=1185876217~%")
(format #t "  for a discussion in which experiments other than~%")
(format #t "  EPR-B experiments are mistaken for EPR-B and analyzed~%")
(format #t "  instead. Perhaps the simulation above is an example of~%")
(format #t "  such a ‘Clauser-correlation experiment’.~%")
(format #t "~%")

(format #t "  ------------------------------------------------------~%")

(format #t "~%")
(format #t "  Here are numbers from an incorrectly done correlation~%")
(format #t "  coefficient calculation—~%")
(format #t "~%")

;;; See
;;; https://en.wikipedia.org/w/index.php?title=CHSH_inequality&oldid=1185876217
(define S-computed 0) ;; For computing a CHSH contrast.
(set! i 1)            ;; For computing a CHSH contrast.

(do ((test-angles bell-test-angles (cdr test-angles)))
    ((null? test-angles))
  (let* ((φ₁ (caar test-angles))
         (φ₂ (cadar test-angles))
         (probs-list (detection-probabilities φ₁ φ₂))
         (computed-correlation (compute-correlation probs-list)))
    (set! S-computed
      ((if (= i 2) - +) S-computed computed-correlation))
    (set! i (+ i 1))
    (format #t "  test angles:  φ₁ = ~A   φ₂ = ~A~%"
            (radians->string φ₁) (radians->string φ₂))
    (do ((patterns pattern-list (cdr patterns)))
        ((null? patterns))
      (let* ((patt (car patterns))
             (prob (cadr (assoc patt probs-list))))
        (format #t "  ~A freq~14,5@F~%" patt (inexact prob))))
    (format #t "     correlation~14,5@F~%" computed-correlation)
    (format #t "~%")))
(format #t "    CHSH S value~14,5@F~%" S-computed)
(format #t "~%")
