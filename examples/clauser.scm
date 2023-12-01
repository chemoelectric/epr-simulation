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

(define θH 0)
(define θV π/2)

(define (detection-probabilities φ₁ φ₂)
  (define pbs₁ (make-pbs φ₁))
  (define pbs₂ (make-pbs φ₂))
  (let-values (((PH₁V₂ PV₁H₂) (photon-pair-probabilities))
               ((PH₁+ PH₁-) (pbs-probabilities pbs₁ θH))
               ((PV₁+ PV₁-) (pbs-probabilities pbs₁ θV))
               ((PH₂+ PH₂-) (pbs-probabilities pbs₂ θH))
               ((PV₂+ PV₂-) (pbs-probabilities pbs₂ θV)))
    `(((H V + +) ,(* PH₁V₂ PH₁+ PV₂+))
      ((H V + -) ,(* PH₁V₂ PH₁+ PV₂-))
      ((H V - +) ,(* PH₁V₂ PH₁- PV₂+))
      ((H V - -) ,(* PH₁V₂ PH₁- PV₂-))
      ((V H + +) ,(* PV₁H₂ PV₁+ PH₂+))
      ((V H + -) ,(* PV₁H₂ PV₁+ PH₂-))
      ((V H - +) ,(* PV₁H₂ PV₁- PH₂+))
      ((V H - -) ,(* PV₁H₂ PV₁- PH₂-)))))

;; (define (correlation-coefficient probabilities)
;;   ;;
;;   ;; The correlation coefficient, assigning +1 to (+) detections and
;;   ;; -1 to (-) detections. With this assignment, the correlation
;;   ;; coefficient equals the covariance.
;;   ;;
;;   (match probabilities
;;     (((_ . P11+) (_ . P11-) (_ . P12+) (_ . P12-)
;;       (_ . P21+) (_ . P21-) (_ . P22+) (_ . P22-))
;;      (let ((P1122++ (* 1/2 P11+ P22+))
;;            (P1122+- (* 1/2 P11+ P22-))
;;            (P1122-+ (* 1/2 P11- P22+))
;;            (P1122-- (* 1/2 P11- P22-))
;;            (P1221++ (* 1/2 P12+ P21+))
;;            (P1221+- (* 1/2 P12+ P21-))
;;            (P1221-+ (* 1/2 P12- P21+))
;;            (P1121-- (* 1/2 P12- P21-)))
;;        (let ((sum-alike (+ P1122++ P1122-- P1221++ P1121--))
;;              (sum-unalike (+ P1122+- P1122-+ P1221+- P1221-+)))
;;          (- sum-alike sum-unalike))))))

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
;; (write (quantum-mechanical-correlation (quantum-mechanical-probabilities π/4 π/8)))(newline)
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
