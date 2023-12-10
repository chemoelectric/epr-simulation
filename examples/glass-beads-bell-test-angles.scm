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
;;; Date first completed : 5 December 2023
;;;
;;; Simulation of a two-channel Bell test with devices that sort glass
;;; beads. The devices behave a lot like Stern-Gerlach magnets. See,
;;; for instance,
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
;;; All the activity is ‘local realistic’. (Surely there is no such
;;; thing as ‘instantaneous action at a distance’.)
;;;

(import (scheme base)
        (scheme write)
        (only (srfi 27) random-real)
        (epr-simulation))

(cond-expand
  (chicken (import (format))
           (import (matchable))))

(define bell-test-angles
  ;; These angles are twice those for an optical experiment.
  `((0 ,π/4)
    (0 ,π3/4)
    (,π/2 ,π/4)
    (,π/2 ,π3/4)))

;; The particles are glass beads. They come in ‘white’ and ‘black’,
;; represented by the symbols 'W and 'B. White beads follow the cos²
;; rule, black beads follow the sin² rule.
(define (bead-pair-source)
  (if (< (random-real) 1/2) (values 'W 'B) (values 'B 'W)))

(define (detection-probabilities φ₁ φ₂)
  ;;
  ;; The probabilities of detections are computed below without
  ;; quantum mechanics. This is the solution demanded by probability
  ;; theory. It is merely the division of possible events into their
  ;; proper proportions, as dictated by the experimental layout.
  ;;
  ;; (There really is NO PHYSICS in the solution of the Bell test. It
  ;; is probability theory, purely. The actual objects and devices can
  ;; be substituted freely, as long as they behave similarly, which
  ;; many do. This implies a quantum mechanical solution is not ‘doing
  ;; physics’, but doing the mathematics by a different method.)
  ;;
  (define splitter₁ (make-splitter 'stern-gerlach φ₁ '(W) '(B)))
  (define splitter₂ (make-splitter 'stern-gerlach φ₂ '(W) '(B)))
  (let-values (((PW₁+ PW₁-) (splitter-probabilities splitter₁ 'W))
               ((PB₁+ PB₁-) (splitter-probabilities splitter₁ 'B))
               ((PW₂+ PW₂-) (splitter-probabilities splitter₂ 'W))
               ((PB₂+ PB₂-) (splitter-probabilities splitter₂ 'B)))
    ;; (W + B +) -- white bead (+) at splitter₁  black bead (+) at splitter₂
    ;; (W + B -) -- white bead (+) at splitter₁  black bead (-) at splitter₂
    ;; etc.
    (let ((probs `(((W + B +) ,(* 1/2 PW₁+ PB₂+))
                   ((W + B -) ,(* 1/2 PW₁+ PB₂-))
                   ((W - B +) ,(* 1/2 PW₁- PB₂+))
                   ((W - B -) ,(* 1/2 PW₁- PB₂-))
                   ((B + W +) ,(* 1/2 PB₁+ PW₂+))
                   ((B + W -) ,(* 1/2 PB₁+ PW₂-))
                   ((B - W +) ,(* 1/2 PB₁- PW₂+))
                   ((B - W -) ,(* 1/2 PB₁- PW₂-)))))

      ;; Sanity check: verify that the probabilities add up to one.
      (check-probabilities (map cadr probs))

      probs)))

(define (simulate-one-event splitter₁ splitter₂)
  (let*-values (((bead₁ bead₂) (bead-pair-source))
                ((detect₁+ _detect₁-) (splitter-activity splitter₁ bead₁))
                ((detect₂+ _detect₂-) (splitter-activity splitter₂ bead₂)))
    ;; (W + B +) -- white bead (+) at splitter₁  black bead (+) at splitter₂
    ;; (W + B -) -- white bead (+) at splitter₁  black bead (-) at splitter₂
    ;; etc.
    `(,bead₁ ,(if detect₁+ '+ '-) ,bead₂ ,(if detect₂+ '+ '-))))

(define *events-per-test-angle* (make-parameter 1000000))

(define (detection-frequencies φ₁ φ₂)
  ;; Simulate events and compute frequencies of the different
  ;; detection patterns.
  (define splitter₁ (make-splitter 'stern-gerlach φ₁ '(W) '(B)))
  (define splitter₂ (make-splitter 'stern-gerlach φ₂ '(W) '(B)))
  (define N (*events-per-test-angle*))
  (define NW+B+ 0) (define NW+B- 0)
  (define NW-B+ 0) (define NW-B- 0)
  (define NB+W+ 0) (define NB+W- 0)
  (define NB-W+ 0) (define NB-W- 0)
  (do ((i 0 (+ i 1)))
      ((= i N))
    (match (simulate-one-event splitter₁ splitter₂)
      ('(W + B +) (set! NW+B+ (+ NW+B+ 1)))
      ('(W + B -) (set! NW+B- (+ NW+B- 1)))
      ('(W - B +) (set! NW-B+ (+ NW-B+ 1)))
      ('(W - B -) (set! NW-B- (+ NW-B- 1)))
      ('(B + W +) (set! NB+W+ (+ NB+W+ 1)))
      ('(B + W -) (set! NB+W- (+ NB+W- 1)))
      ('(B - W +) (set! NB-W+ (+ NB-W+ 1)))
      ('(B - W -) (set! NB-W- (+ NB-W- 1)))))
  ;; (W + B +) -- white bead (+) at splitter₁  black bead (+) at splitter₂
  ;; (W + B -) -- white bead (+) at splitter₁  black bead (-) at splitter₂
  ;; etc.
  `(((W + B +) ,(/ NW+B+ N))
    ((W + B -) ,(/ NW+B- N))
    ((W - B +) ,(/ NW-B+ N))
    ((W - B -) ,(/ NW-B- N))
    ((B + W +) ,(/ NB+W+ N))
    ((B + W -) ,(/ NB+W- N))
    ((B - W +) ,(/ NB-W+ N))
    ((B - W -) ,(/ NB-W- N))))

(define (estimate-correlation φ₁ φ₂ detection-freqs)
  ;; Use detection frequencies to estimate the value of -cos(φ₁-φ₂).
  (define (get-freq pattern)
    (cadr (assoc pattern detection-freqs)))
  (estimate-pair-correlation
   φ₁ φ₂ 'stern-gerlach 'complementary
   `(,(get-freq '(W + B +)) ,(get-freq '(W + B -))
     ,(get-freq '(W - B +)) ,(get-freq '(W - B -))
     ,(get-freq '(B + W +)) ,(get-freq '(B + W -))
     ,(get-freq '(B - W +)) ,(get-freq '(B - W -)))))

(define pattern-list
  '((W + B +) (W + B -) (W - B +) (W - B -)
    (B + W +) (B + W -) (B - W +) (B - W -)))

(format #t "~%")
(format #t "  A simulation of glass beads passing through two-channel~%")
(format #t "  devices that behave similarly to Stern-Gerlach magnets.~%")
(format #t "~%")
(format #t "  legend:~%")
(format #t "    (W + B +)  white bead in (+) channel of splitter₁,~%")
(format #t "               black bead in (+) channel of splitter₂,~%")
(format #t "    (W + B -)  white bead in (+) channel of splitter₁,~%")
(format #t "               black bead in (-) channel of splitter₂, etc.~%")
(format #t "~%")
(format #t "  White beads follow the cos² rule.~%")
(format #t "  Black beads follow the sin² rule.~%")
(format #t "  Nominal correlation is −cos(φ₁ − φ₂).~%")

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
         (nominal-correlation (- (cos (- φ₁ φ₂))))
         (estimated-correlation
          (estimate-correlation φ₁ φ₂ freqs-list)))
    (set! S-nominal
      ((if (= i 2) - +) S-nominal nominal-correlation))
    (set! S-estimated
      ((if (= i 2) - +) S-estimated estimated-correlation))
    (set! i (+ i 1))
    (format #t "  test angles:  φ₁ = ~A   φ₂ = ~A~%"
            (radians->string φ₁) (radians->string φ₂))
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
(format #t "  See https://en.wikipedia.org/w/index.php?title=CHSH_inequality&oldid=1185876217~%")
(format #t "  for a page making the false claim this result cannot be~%")
(format #t "  gotten without ‘quantum entanglement’. Note ‘quantum~%")
(format #t "  correlation’ is simply the correlation coefficient when~%")
(format #t "  computed by quantum mechanics rather than other means.~%")
(format #t "~%")
(format #t "  Of course, ‘quantum entanglement’ of glass beads would~%")
(format #t "  be laughed out of the room. However, the word problem~%")
(format #t "  was always equivalent to this one about glass beads!~%")
(format #t "~%")
