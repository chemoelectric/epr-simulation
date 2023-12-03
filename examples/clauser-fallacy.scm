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
                ;; phot₁ simply gets counted. (The following is
                ;; equivalent to a PBS set to zero.)
                ((detect₁+) (eq? (photon->symbol phot₁) 'H))
                ;; But phot₂ goes through two PBSes.
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

  ;; The first PBS outputs photons whose polarization angles we need
  ;; to take into account.
  (define pbs₁ (make-pbs φ₁ φ₁ (+ φ₁ π/2)))

  ;; The second PBS outputs photons into photodetectors. We can ignore
  ;; polarization.
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
         (nominal-correlation (compute-correlation probs-list))
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
    (format #t "   \"correlation\"~14,5@F~14,5@F~%"
            nominal-correlation estimated-correlation)
    (format #t "~%")))
(format #t "                       nominal     simulated~%")
(format #t "    CHSH S value~14,5@F~14,5@F~%"
        S-nominal S-estimated)
(format #t "~%")
(format #t "  (See https://en.wikipedia.org/w/index.php?title=CHSH_inequality&oldid=1185876217~%")
(format #t "  for a discussion in which an experiment of the kind~%")
(format #t "  simulated here is mistaken for an EPR-B Bell test.~%")
(format #t "  Notice the simulation does NOT produce the predicted~%")
(format #t "  frequencies of events, but DOES approximate the~%")
(format #t "  predicted correlations. This happens because Clauser-~%")
(format #t "  style calculations of the correlation coefficient are~%")
(format #t "  done incorrectly, and give answers for an experiment~%")
(format #t "  such as we simulate here, instead of the Bell test.)~%")
(format #t "~%")
