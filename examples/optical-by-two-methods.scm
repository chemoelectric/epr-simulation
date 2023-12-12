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
;;; Date first completed : ????????????????????????????
;;;
;;; Simulation of a two-channel optical Bell test by two different
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

;; (define (photon->symbol phot)
;;   (if (< (photon-polarization-angle phot) 0.0001) 'H 'V))

;; (define (simulate-one-event pbs₁ pbs₂)
;;   (let*-values (((phot₁ phot₂) (photon-pair-source θH θV))
;;                 ((detect₁+ _detect₁-) (pbs-activity pbs₁ phot₁))
;;                 ((detect₂+ _detect₂-) (pbs-activity pbs₂ phot₂)))
;;     ;; (H + V +) -- horiz (+) at pbs₁  vert (+) at pbs₂
;;     ;; (H + V -) -- horiz (+) at pbs₁  vert (-) at pbs₂
;;     ;; etc.
;;     `(,(photon->symbol phot₁) ,(if detect₁+ '+ '-)
;;       ,(photon->symbol phot₂) ,(if detect₂+ '+ '-))))

;; (define *events-per-test-angle* (make-parameter 1000000))

;; (define (detection-frequencies φ₁ φ₂)
;;   ;; Simulate events and compute frequencies of the different
;;   ;; detection patterns.
;;   (define pbs₁ (make-pbs φ₁))
;;   (define pbs₂ (make-pbs φ₂))
;;   (define N (*events-per-test-angle*))
;;   (define NH+V+ 0) (define NH+V- 0)
;;   (define NH-V+ 0) (define NH-V- 0)
;;   (define NV+H+ 0) (define NV+H- 0)
;;   (define NV-H+ 0) (define NV-H- 0)
;;   (do ((i 0 (+ i 1)))
;;       ((= i N))
;;     (match (simulate-one-event pbs₁ pbs₂)
;;       ('(H + V +) (set! NH+V+ (+ NH+V+ 1)))
;;       ('(H + V -) (set! NH+V- (+ NH+V- 1)))
;;       ('(H - V +) (set! NH-V+ (+ NH-V+ 1)))
;;       ('(H - V -) (set! NH-V- (+ NH-V- 1)))
;;       ('(V + H +) (set! NV+H+ (+ NV+H+ 1)))
;;       ('(V + H -) (set! NV+H- (+ NV+H- 1)))
;;       ('(V - H +) (set! NV-H+ (+ NV-H+ 1)))
;;       ('(V - H -) (set! NV-H- (+ NV-H- 1)))))
;;   ;; (H + V +) -- horiz (+) at pbs₁  vert (+) at pbs₂
;;   ;; (H + V -) -- horiz (+) at pbs₁  vert (-) at pbs₂
;;   ;; etc.
;;   `(((H + V +) ,(/ NH+V+ N))
;;     ((H + V -) ,(/ NH+V- N))
;;     ((H - V +) ,(/ NH-V+ N))
;;     ((H - V -) ,(/ NH-V- N))
;;     ((V + H +) ,(/ NV+H+ N))
;;     ((V + H -) ,(/ NV+H- N))
;;     ((V - H +) ,(/ NV-H+ N))
;;     ((V - H -) ,(/ NV-H- N))))

;; (define (estimate-correlation φ₁ φ₂ detection-freqs)
;;   ;; Use detection frequencies to estimate the value of -cos(2(φ₁-φ₂).
;;   (define (get-freq pattern)
;;     (cadr (assoc pattern detection-freqs)))
;;   (estimate-pair-correlation
;;    φ₁ φ₂ 'optical 'complementary
;;    `(,(get-freq '(H + V +)) ,(get-freq '(H + V -))
;;      ,(get-freq '(H - V +)) ,(get-freq '(H - V -))
;;      ,(get-freq '(V + H +)) ,(get-freq '(V + H -))
;;      ,(get-freq '(V - H +)) ,(get-freq '(V - H -)))))

;; (define pattern-list
;;   '((H + V +) (H + V -) (H - V +) (H - V -)
;;     (V + H +) (V + H -) (V - H +) (V - H -)))

;; (format #t "~%")
;; (format #t "  Simulation of a two-channel optical Bell test.~%")
;; (format #t "~%")
;; (format #t "  legend:~%")
;; (format #t "    (H + V +)  horizontal photon in (+) channel of pbs₁,~%")
;; (format #t "               vertical photon in (+) channel of pbs₂,~%")
;; (format #t "    (H + V -)  horizontal photon in (+) channel of pbs₁,~%")
;; (format #t "               vertical photon in (-) channel of pbs₂, etc.~%")

;; ;;; See
;; ;;; https://en.wikipedia.org/w/index.php?title=CHSH_inequality&oldid=1185876217
;; ;;; Be aware that the contents of the page is mostly wrong. That is
;; ;;; proven by this simulation. This contrast calculation is included
;; ;;; to emphasize our proof.
;; (define S-nominal 0)   ;; For computing a CHSH contrast.
;; (define S-estimated 0) ;; For computing a CHSH contrast.
;; (define i 1)           ;; For computing a CHSH contrast.

;; (format #t "~%")
;; (do ((test-angles bell-test-angles (cdr test-angles)))
;;     ((null? test-angles))
;;   (let* ((φ₁ (caar test-angles))
;;          (φ₂ (cadar test-angles))
;;          (probs-list (detection-probabilities φ₁ φ₂))
;;          (freqs-list (detection-frequencies φ₁ φ₂))
;;          (nominal-correlation (- (cos (* 2 (- φ₁ φ₂)))))
;;          (estimated-correlation
;;           (estimate-correlation φ₁ φ₂ freqs-list)))
;;     (set! S-nominal
;;       ((if (= i 2) - +) S-nominal nominal-correlation))
;;     (set! S-estimated
;;       ((if (= i 2) - +) S-estimated estimated-correlation))
;;     (set! i (+ i 1))
;;     (format #t "  test angles:  φ₁ = ~A   φ₂ = ~A~%"
;;             (radians->string φ₁) (radians->string φ₂))
;;     (format #t "                       nominal     simulated~%")
;;     (do ((patterns pattern-list (cdr patterns)))
;;         ((null? patterns))
;;       (let* ((patt (car patterns))
;;              (prob (cadr (assoc patt probs-list)))
;;              (freq (cadr (assoc patt freqs-list))))
;;         (format #t "  ~A freq~14,5@F~14,5@F~%"
;;                 patt (inexact prob) (inexact freq))))
;;     (format #t "     correlation~14,5@F~14,5@F~%"
;;             nominal-correlation estimated-correlation)
;;     (format #t "~%")))
;; (format #t "                       nominal     simulated~%")
;; (format #t "    CHSH S value~14,5@F~14,5@F~%"
;;         S-nominal S-estimated)
;; (format #t "~%")
;; (format #t "  (See https://en.wikipedia.org/w/index.php?title=CHSH_inequality&oldid=1185876217~%")
;; (format #t "  for a page making the false claim this result cannot be~%")
;; (format #t "  gotten without ‘quantum entanglement’. Note ‘quantum~%")
;; (format #t "  correlation’ is simply the correlation coefficient when~%")
;; (format #t "  computed by quantum mechanics rather than other means.)~%")
;; (format #t "~%")

(define t1 (tensor-normalize
            '((1.2 . "a,b,c") (2.3 . "b,a,c")  (3.4 . "c,a,b"))))
(write t1)(newline)
(define t1.0 (tensor./ t1 0))
(define t1.1 (tensor./ t1 1))
(define t1.2 (tensor./ t1 2))
(write t1.0)(newline)
(write t1.1)(newline)
(write t1.2)(newline)
