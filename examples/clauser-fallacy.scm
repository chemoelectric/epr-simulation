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
;;; Date first completed : ??????????????????????
;;;
;;; Simulation of what the mistaken analysis of John Clauser and
;;; others ACTUALLY seems to represent. See, for instance,
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
;;; Bell’s original notation) is NOT an expectation for an EPR-B
;;; experiment.
;;;
;;;    IT IS AN EXPECTATION FOR A SINGLE PHOTON PASSING THROUGH TWO
;;;    POLARIZING BEAM SPLITTERS IN SERIES!
;;;
;;; This should be intuitive once pointed out. It explains why, in the
;;; analyses of Clauser, the correlations are supposed to disappear
;;; when a Bell test angle is changed from zero to π/4. An ideal
;;; polarizing beam splitter set to zero could literally be left out
;;; of the apparatus, without effect. One will get the same results as
;;; in the EPR-B experiment. A polarizing beam splitter set to π/4, on
;;; the other hand, will obstruct the particle beam.
;;;

(import (scheme base)
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

      ;; Sanity check: verify that the probabilities add up to 1.
      (let ((sum (apply + (map cadr probs))))
        (unless (<= (abs (- sum 1)) (* 500 fl-epsilon))
          (error "detection-probabilities: probabilities do not add up to 1")))

      probs)))

(define (compute-correlation probabilities)
  ;;
  ;; The correlation coefficient, assigning +1 to (+) detections and
  ;; -1 to (-) detections. With this assignment, the correlation
  ;; coefficient equals the covariance.
  ;;
  ;; This will be a ‘Clauser-style’ calculation, for polarizing beam
  ;; splitters in series.
  ;;
  ;; The orthodoxy in physics has mistaken this as a calculation for
  ;; the Bell test. The correct calculation for the Bell test,
  ;; however, obviously has to come out the same as what quantum
  ;; mechanics calculates. This is required by the principle that all
  ;; mathematical methods must produce the same results.
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

                ((detect₁+ _detect₁-)
                 ;; phot₁ simply gets counted. (The following is
                 ;; equivalent to a PBS set to zero.)
                 (if (eq? (photon->symbol phot₁) 'H)
                     (values #t #f)
                     (values #f #t)))

                ;; But phot₂ goes through two PBSes.
                ((phot₂+ phot₂-) (pbs-activity pbs₁ phot₂))
                ((detect₂++ _detect₂+-) (pbs-activity pbs₂ phot₂+))
                ((detect₂-+ _detect₂--) (pbs-activity pbs₂ phot₂-)))

    (let ((detect₂+ (or detect₂++ detect₂-+))
          (detect₂- (or detect₂+- detect₂--)))

      ;; (H + V +) -- horiz (+) at pbs₁  vert (+) at pbs₂
      ;; (H + V -) -- horiz (+) at pbs₁  vert (-) at pbs₂
      ;; etc.
      `(,(photon->symbol phot₁) ,(if detect₁+ '+ '-)
        ,(photon->symbol phot₂) ,(if detect₂+ '+ '-)))))


;; (define (simulate-one-event φ₁ φ₂)
;;   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FIXME: Simulate the Clauser experiment.
;;   (define pbs₁ (make-pbs φ₁))
;;   (define pbs₂ (make-pbs φ₂))
;;   (let*-values (((phot₁ phot₂) (photon-pair-source θ₁ θ₂))
;;                 ((detect₁+ _detect₁-) (pbs-activity pbs₁ phot₁))
;;                 ((detect₂+ _detect₂-) (pbs-activity pbs₂ phot₂)))
;;     `((,phot₁ ,(if detect₁+ '+ '-))
;;       (,phot₂ ,(if detect₂+ '+ '-)))))

;; (write (quantum-mechanical-probabilities 0 π/8))(newline)
(write (compute-correlation (detection-probabilities 0 π/8)))(newline)
(write (compute-correlation (detection-probabilities π/4 π/8)))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
;; (write (simulate-one-event 0 π/8))(newline)
