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
;;; Date first completed : 19 December 2023
;;;
;;; Analysis of a two-channel optical Bell test by two different
;;; methods.
;;; https://en.wikipedia.org/w/index.php?title=CHSH_inequality&oldid=1185876217
;;;
;;; The analysis there is wrong, and the experiments run by Aspect et
;;; alia are meaningless.
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
;;; Here the probabilities of outcomes of the Bell test are calculated
;;; both by probability theory and by tensors (quantum mechanics).
;;;

(import (scheme base)
        (scheme write)
        (srfi 1)                        ; R⁷RS-large (scheme list)
        (epr-simulation))

(cond-expand
  (chicken (import (format))
           (import (matchable))))

(define my-test-angles
  `((0 ,π/8)
    (0 ,π/4)
    (0 ,π3/8)
    (0 ,π/2)
    (,π/4 ,0)
    (,π/4 ,π/8)
    (,π/4 ,π3/8)
    (,π/4 ,π/2)))

(define θH 0)
(define θV π/2)

(define (detection-probabilities-by-probability-theory φ₁ φ₂)
  ;;
  ;; The probabilities of detections are computed below without
  ;; quantum mechanics. This is the solution demanded by probability
  ;; theory. It is merely the division of possible events into their
  ;; proper proportions.
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

(define (system-state φ₁ φ₂)
  ;;
  ;; Return the tensor for the system state according to quantum
  ;; mechanics. Although some attribute physical meaning to this
  ;; expression, it is actually an obfuscated form of the probability
  ;; calculation. Points in the linear space represent states of a
  ;; computation, not of physical entities.
  ;;
  ;; Because quantum mechanics is a method of calculation and not a
  ;; physical theory, we can feel free to include whatever information
  ;; is available. Thus the final system state will include the
  ;; original photon polarization angles.
  ;;
  (define pbs₁ (make-pbs φ₁))
  (define pbs₂ (make-pbs φ₂))
  (let-values (((PH₁V₂ PV₁H₂) (photon-pair-probabilities)))
    (let ((pbs₁H (tensor./ (pbs-tensor pbs₁ θH) '(0 1)))
          (pbs₁V (tensor./ (pbs-tensor pbs₁ θV) '(0 1)))
          (pbs₂H (tensor./ (pbs-tensor pbs₂ θH) '(0 1)))
          (pbs₂V (tensor./ (pbs-tensor pbs₂ θV) '(0 1))))
      (tensor+ (tensor* (sqrt PH₁V₂) (tensor.* pbs₁H pbs₂V))
               (tensor* (sqrt PV₁H₂) (tensor.* pbs₁V pbs₂H))))))

(define (detection-probabilities-by-quantum-mechanics φ₁ φ₂)
  (map (lambda (term) `(,(tensor-basis->HV-symbols (cdr term))
                        ,(car term)))
       (tensor->probabilities (system-state φ₁ φ₂))))

(define (angle->symbol angle)
  (if (< angle 0.0001) 'H 'V))

(define (tensor-basis->HV-symbols basis)
  `(,(angle->symbol (string->radians (tensor-basis-ref basis 0)))
    ,(string->symbol (tensor-basis-ref basis 1))
    ,(angle->symbol (string->radians (tensor-basis-ref basis 2)))
    ,(string->symbol (tensor-basis-ref basis 3))))

(define pattern-list
  '((H + V +) (H + V -) (H - V +) (H - V -)
    (V + H +) (V + H -) (V - H +) (V - H -)))

(format #t "~%")
(format #t "  Calculation of probabilities of outcomes of a~%")
(format #t "  two-channel optical Bell test by two methods:~%")
(format #t "  ordinary probability theory and quantum mechanics.~%")
(format #t "~%")
(format #t "  legend:~%")
(format #t "    (H + V +)  horizontal photon in (+) channel of pbs₁,~%")
(format #t "               vertical photon in (+) channel of pbs₂,~%")
(format #t "    (H + V -)  horizontal photon in (+) channel of pbs₁,~%")
(format #t "               vertical photon in (-) channel of pbs₂, etc.~%")

(do ((test-angles my-test-angles (cdr test-angles)))
    ((null? test-angles))
  (let* ((φ₁ (caar test-angles))
         (φ₂ (cadar test-angles))
         (pprobs
          (detection-probabilities-by-probability-theory φ₁ φ₂))
         (qprobs
          (detection-probabilities-by-quantum-mechanics φ₁ φ₂)))
    (format #t "~%")
    (format #t "  test angles:  φ₁ = ~A   φ₂ = ~A~%"
            (radians->string φ₁) (radians->string φ₂))
    (format #t "              probability       quantum~%")
    (format #t "                   theory     mechanics~%")
    (do ((patterns pattern-list (cdr patterns)))
        ((null? patterns))
      (let* ((patt (car patterns))
             (pprob (cadr (assoc patt pprobs)))
             (qprob (cadr (assoc patt qprobs))))
        (format #t "    ~A~12,5@F~14,5@F~%"
                patt (inexact pprob) (inexact qprob))))))
(format #t "~%")
