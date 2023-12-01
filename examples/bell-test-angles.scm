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
;;; All the activity is ‘local realistic’ (as opposed to non-local
;;; preposterous). No educated person should give ‘action at a
;;; distance’ even a moment’s credence. It is logically inconsistent
;;; with the bulk of human knowledge.
;;;

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

(define θH 0)
(define θV π/2)

(define (detection-probabilities φ₁ φ₂)
  ;;
  ;; The probabilities of detections is solved below easily, without
  ;; quantum mechanics.
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
  ;; amplitudes. There is no physical meaning. It is merely a way to
  ;; do probability calculations, employed by people who do not know
  ;; probability theory.
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


;; (define (quantum-mechanical-correlation probabilities)
;;   ;;
;;   ;; The correlation coefficient, assigning +1 to (+) detections and
;;   ;; -1 to (-) detections. With this assignment, the correlation
;;   ;; coefficient equals the covariance.
;;   ;;
;;   ;; By the way, one could also use the closed solution
;;   ;; -cos(2(θ₁-θ₂))=-(cos²(θ₁-θ₂)-sin²(θ₁-θ₂)). This solution has been
;;   ;; derived from quantum mechanics—but ALSO by the author, using ONLY
;;   ;; probability theory, WITHOUT quantum mechanics. The closed
;;   ;; solution ALSO follows from classical coherence theory, if one
;;   ;; assumes any of various statistical mechanical points of view.
;;   ;;
;;   ;; Orthodox ‘quantum’ physics is astoundingly WRONG about the Bell
;;   ;; test, despite having awarded Nobel Prizes for it. You are
;;   ;; witnessing that wrongness here.
;;   ;;
;;   (match probabilities
;;     (((_ . P11+) (_ . P11-) (_ . P12+) (_ . P12-)
;;       (_ . P21+) (_ . P21-) (_ . P22+) (_ . P22-))
;;      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FIXME: THIS IS WRONG!!!!!!!!!!!!!!!!!!!!!!!!! IT IS THE CLAUSER MISTAKE!!!!!!!!!
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

(write (detection-probabilities 0 π/8))(newline)
(write (map (lambda (x) (cons (car x) (inexact (cadr x)))) (detection-frequencies 0 π/8)))(newline)
;; (write (quantum-mechanical-correlation (quantum-mechanical-probabilities π/4 π/8)))(newline)
(write (simulate-one-event (make-pbs 0) (make-pbs π/8)))(newline)
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
