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

(format #t "~%")
(do ((test-angles bell-test-angles (cdr test-angles)))
    ((null? test-angles))
  (let* ((φ₁ (caar test-angles))
         (φ₂ (cadar test-angles))
         (probs-list (detection-probabilities φ₁ φ₂))
         (freqs-list (detection-frequencies φ₁ φ₂)))
    (format #t "~2Ttest angles:  ~16Tφ₁ = ~A~33Tφ₂ = ~A~%"
            (angle->string φ₁) (angle->string φ₂))
    (format #t "~4Tdetection~18Tnominal~32Tsimulated~%")
    (format #t "~4Tpattern~18Tprobability~32Tfrequency~%")
    (do ((patterns pattern-list (cdr patterns)))
        ((null? patterns))
      (let* ((patt (car patterns))
             (prob (cadr (assoc patt probs-list)))
             (freq (cadr (assoc patt freqs-list))))
        (format #t "~4T~A~18T~8,5F~32T~8,5F~%"
                patt (inexact prob) (inexact freq))))
    (format #t "~%")))
