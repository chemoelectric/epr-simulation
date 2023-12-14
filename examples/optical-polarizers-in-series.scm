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
;;; Date first completed : 14 December 2023
;;;
;;;
;;; The following simulation is of one possible arrangement of two
;;; polarizers in series, and it gives a CHSH ‘contrast’ |S|=√2 < 2.
;;;
;;; This simulation perhaps is suggestive of mistaken calculations
;;; that underlie such fallacies as the CHSH inequality.
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
(format #t "  Be forewarned the page is full of erroneous analysis.~%")
(format #t "~%")
