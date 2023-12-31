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

(define-library (epr-simulation)

  (export π/180 π/8 π/4 π3/8 π/2 π3/4 π)

  (export degrees->radians radians->degrees)
  (export radians->string string->radians) ; Using the prefix π or π×.

  (export approx=)

  ;; Support for vectors and tensors (ordered tuples of vectors,
  ;; forming vectors themselves) as lists of cons-pairs. Each car is a
  ;; number representing the amplitude and each cdr is a string
  ;; representing the basis tensor. This is a simplified notation for
  ;; all that is needed to handle EPR-B by quantum mechanics.
  (export tensor+                  ; Add tensors or collect terms.
          tensor*                  ; Scalar multiplication.
          tensor.*                 ; Lengthen the tuples.
          tensor./                 ; Shorten or rearrange the tuples.
          tensor->probabilities    ; Convert from amplitudes to
                                        ; probabilities.
          tensor->amplitudes       ; Convert from probabilities to
                                        ; amplitudes.
          tensor-normalize         ; Normalize a tensor to probability
                                        ; one. Use with care.
          tensor-basis-split       ; Split a basis tensor string.
          tensor-basis-ref)        ; Return the i’th part of a split.

  (export <photon>
          ;; A plane-polarized photon.
          make-photon
          photon?
          photon-polarization-angle set-photon-polarization-angle!
          make-photon-polarization-angle-changer
          *photon-pair-probability* ; Parameter, default 0.5
          photon-pair-probabilities ; Two values, default 0.5, 0.5
          photon-pair-source)       ; Two photons, complementary pair.

  (export <pbs>
          ;; Polarizing beam-splitter with plane-polarized output.
          make-pbs
          pbs?
          pbs-angle          ; Angle wrt the incident photon.
          pbs-probabilities  ; Probabilities given an incident photon.
          pbs-tensor         ; Quantum mechanical representation.
          pbs-activity)      ; Activity given an incident photon.

  (export <splitter>
          ;; Any of many possible two-channel devices that sort
          ;; objects into channels according to cos² and sin² rules.
          make-splitter
          splitter?
          splitter-behavior           ; 'optical or 'stern-gerlach
          splitter-angle              ; Angle wrt the incident object.
          splitter-cos²-rule          ; Objects that obey cos² rule.
          splitter-cos2-rule          ; Synonym.
          splitter-sin²-rule          ; Particles that obey sin² rule.
          splitter-sin2-rule          ; Synonym.
          splitter-probabilities      ; Probabilities, given an
                                        ; incident object.
          splitter-activity)          ; Activity, given an incident
                                        ; object.

  (export estimate-pair-correlation)  ; Estimate correlation from
                                        ; detection frequencies.

  (export check-probabilities)        ; Check if a list of
                                        ; probabilities adds up
                                        ; approximately to one.

  (import (scheme base)
          (scheme case-lambda)
          (scheme complex)
          (scheme cxr)
          (scheme inexact)
          (srfi 1)                      ; R⁷RS-large (scheme list)
          (only (srfi 27) random-real)
          (only (srfi 132) list-sort)   ; R⁷RS-large (scheme sort)
          (only (srfi 144) fl-epsilon)) ; R⁷RS-large (scheme flonum)

  ;; It is necessary to use a string library that understands some
  ;; Unicode. For CHICKEN Scheme, this will suffice.
  (import (srfi 152))                   ; String library (reduced).

  (cond-expand
    (chicken (import (only (chicken base) define-record-printer)
                     (scheme write)))
    (else))

  (begin

    (cond-expand
      (chicken (define (write-angle angle port)
                 (display (radians->string angle) port)))
      (else))

    ;; OEIS A019685
    (define π/180 0.0174532925199432957692369076848861271344287188854172545609719144017100911460344944368224156963450948)
    
    ;; OEIS A019675
    (define π/8 0.392699081698724154807830422909937860524646174921888227621868074038477050785776124828504353167764633)

    ;; OEIS A003881
    (define π/4 0.785398163397448309615660845819875721049292349843776455243736148076954101571552249657008706335529266995537)

    ;; OEIS A093828
    (define π3/8 1.17809724509617246442349126872981358157393852476566468286560422211543115235732837448551305950329390049)

    ;; OEIS A019669
    (define π/2 1.57079632679489661923132169163975144209858469968755291048747229615390820314310449931401741267105853)

    ;; OEIS A177870
    (define π3/4 2.356194490192344928846982537459627163147877049531329365731208444230862304714656748971026119006587800986611)

    ;; OEIS A000796
    (define π 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214)

    (define (degrees->radians x) (* x π/180))
    (define (radians->degrees x) (/ x π/180))

    (define radians->string
      (case-lambda
        ((angle) (radians->string angle 1000))
        ((angle
          tolerance) ;; In multiples of one half of fl-epsilon.
         (let* ((angle/π (/ angle π))
                (angle/π*64 (* 64 angle/π))
                (iangle/π*64 (round angle/π*64))
                (diff (abs (- angle/π*64 iangle/π*64)))
                (exact-enough (<= (+ diff diff)
                                  (* tolerance fl-epsilon
                                     (abs angle/π*64)))))
           (string-append
            "π×" (let ((angle/π
                        (if exact-enough (exact angle/π) angle/π)))
                   (number->string angle/π)))))))

    (define length-of-π× (string-length "π×"))
    (define length-of-π (string-length "π"))

    (define (string->radians s)
      (cond ((string-prefix? "π×" s)
             (* π (string->number (string-drop s length-of-π×))))
            ((string-prefix? "π" s)
             (* π (string->number (string-drop s length-of-π))))
            (else (string->number s))))

    (define approx=
      (case-lambda
        ((x y) (approx= x y 1000))
        ((x y tolerance) ;; In multiples of one half of fl-epsilon.
         (unless (real? x)
           (error "approx=: expected a real number" x))
         (unless (real? y)
           (error "approx=: expected a real number" y))
         (unless (positive? tolerance)
           (error "approx=: expected a positive number" tolerance))
         (let* ((tol (* tolerance fl-epsilon))
                (tol (* tol (max (abs x) (abs y) (/ tol 2)))))
           (<= (* 2 (abs (- x y))) tol)))))

    (define (%%check-term caller term)
      (unless (and (pair? term)
                   (number? (car term))
                   (string? (cdr term)))
        (error (append caller ": expected a (number . string) pair")
               term)))

    (define (%%check-tensor caller tensor)
      (unless (list? tensor)
        (error (append caller ": expected a list") tensor))
      (for-each (lambda (term) (%%check-term caller term)) tensor))

    (define (%%combine-terms tensor)
      (let ((tensor (list-sort (lambda (term1 term2)
                                 (string>=? (cdr term1) (cdr term2)))
                               tensor)))
        (let loop ((p tensor)
                   (q '()))
          (cond ((null? p) q)
                ((null? q) (loop (cdr p) (cons (car p) q)))
                ((string=? (cdar p) (cdar q))
                 (loop (cdr p) (cons (cons (+ (caar p) (caar q))
                                           (cdar p))
                                     (cdr q))))
                (else (loop (cdr p) (cons (car p) q)))))))

    (define (%%add-amplitude-as-probability term sum)
      (+ sum (square (magnitude (car term)))))

    (define (tensor-normalize tensor)
      ;; Normalize a tensor so its amplitudes represent a total
      ;; probability of one.
      (let* ((tensor (tensor+ tensor))
             (sum (fold %%add-amplitude-as-probability 0 tensor))
             (denom (sqrt sum)))
        (map (lambda (term) (cons (/ (car term) denom) (cdr term)))
             tensor)))

    (define (tensor+ tensor . tensors)
      ;; Add tensors, or combine terms in one tensor.
      (%%check-tensor "tensor+" tensor)
      (for-each (lambda (t) (%%check-tensor "tensor+" t)) tensors)
      (let ((tensor (concatenate (cons tensor tensors))))
        (%%combine-terms tensor)))

    (define (tensor* scalar tensor)
      (if (and (complex? tensor) (pair? scalar))
          (tensor* tensor scalar)
          (begin
            (unless (complex? scalar)
              (error "tensor*: expected a complex number" scalar))
            (%%check-tensor "tensor*" tensor)
            (map (lambda (term)
                   `(,(* scalar (car term)) . ,(cdr term)))
                 tensor))))

    (define (%%tensor.* caller tensor1 tensor2)
      (%%check-tensor caller tensor2)
      (let loop ((t1 tensor1)
                 (t2 tensor2)
                 (p '()))
        (cond ((null? t1) (tensor+ p))
              ((null? t2) (loop (cdr t1) tensor2 p))
              (else (let ((ampl_t1 (caar t1))
                          (basis_t1 (cdar t1))
                          (ampl_t2 (caar t2))
                          (basis_t2 (cdar t2)))
                      (let ((new-term
                             (cons (* ampl_t1 ampl_t2)
                                   (string-append basis_t1 ","
                                                  basis_t2))))
                        (loop t1 (cdr t2) (cons new-term p))))))))

    (define (tensor.* tensor1 tensor2 . tensors)
      (define caller "tensor.*")
      (define (.* t1 t2) (%%tensor.* caller t1 t2))
      (%%check-tensor caller tensor1)
      (let loop ((p tensors)
                 (q (.* tensor1 tensor2)))
        (if (null? p) q (loop (cdr p) (.* q (car p))))))

    (define (tensor./ tensor i-list)
      ;; Extract the (i₁ i₂ ... iₙ)’th vectors (starting from i=0)
      ;; from the ordered tuple. This is done with probabilities
      ;; rather than amplitudes, and returns positive real amplitudes.
      (define (subcomponent v*-string)
        (let* ((v* (tensor-basis-split v*-string))
               (v* (map (lambda (i) (list-ref v* i)) i-list)))
          (string-join v* ",")))
      (%%check-tensor "tensor./" tensor)
      (unless (and (pair? i-list)
                   (list? i-list)
                   (every integer? i-list)
                   (not (any negative? i-list)))
        (error "tensor./: expected list of indices" i-list))
      (let* ((probs
              (map (lambda (term)
                     `(,(square (magnitude (car term)))
                       . ,(subcomponent (cdr term))))
                   tensor))
             (probs (%%combine-terms probs)))
        (map (lambda (term) `(,(sqrt (car term)) . ,(cdr term)))
             probs)))

    (define (tensor->probabilities tensor)
      (map (lambda (term) `(,(square (magnitude (car term)))
                            . ,(cdr term)))
           tensor))

    (define (tensor->amplitudes tensor)
      (map (lambda (term) `(,(sqrt (car term)) . ,(cdr term)))
           tensor))

    (define (tensor-basis-split basis-tensor-string)
      (string-split basis-tensor-string ","))

    (define (tensor-basis-ref basis-tensor-string i)
      (list-ref (tensor-basis-split basis-tensor-string) i))

    (define-record-type <photon>
      ;; A photon is a light pulse of unit intensity. It has some
      ;; polarization angle orthogonal to its direction of motion. The
      ;; polarization angle is given by a real number.
      (%%make-photon angle)
      photon?
      (angle photon-polarization-angle
             %%set-photon-polarization-angle!))

    (define (make-photon angle)
      (unless (real? angle)
        (error "make-photon: expected a real number" angle))
      (%%make-photon angle))

    (define (set-photon-polarization-angle! photon θ)
      (unless (real? θ)
        (error
         "set-photon-polarization-angle!: expected a real number" θ))
      (%%set-photon-polarization-angle! photon θ))

    (define (make-photon-polarization-angle-changer proc)
      (lambda (phot)
        (make-photon (proc (photon-polarization-angle phot)))))

    (cond-expand
      (chicken (define-record-printer (<photon> rectype port)
                 (display "<photon " port)
                 (write-angle (photon-polarization-angle rectype) port)
                 (display ">" port)))
      (else))

    (define *photon-pair-probability*
      (make-parameter
       0.5
       (lambda (p₁)
         (unless (and (real? p₁) (<= 0 p₁) (<= p₁ 1))
           (error
            "*photon-pair-probability*: expected real number in [0,1]"
            p₁))
         p₁)))

    (define (photon-pair-probabilities)
      (let* ((p₁ (*photon-pair-probability*))
             (p₂ (- 1 p₁)))
        (values p₁ p₂)))

    (define (photon-pair-source angle1 angle2)
      (call-with-values photon-pair-probabilities
        (lambda (p₁ _p₂)
          (let* ((lessthan (< (random-real) p₁))
                 (θ₁ (if lessthan angle1 angle2))
                 (θ₂ (if lessthan angle2 angle1)))
            (values (make-photon θ₁) (make-photon θ₂))))))

    (define-record-type <pbs>
      ;; Polarizing beam-splitter with plane-polarized output. (A
      ;; two-channel polarizer.)
      (%%make-pbs angle)
      pbs?
      (angle pbs-angle))

    (define (make-pbs angle)
      (unless (real? angle)
        (error "make-pbs: expected a real number" angle))
      (%%make-pbs angle))

    (cond-expand
      (chicken (define-record-printer (<pbs> rectype port)
                 (display "<pbs " port)
                 (write-angle (pbs-angle rectype) port)
                 (display ">" port)))
      (else))

    (define (pbs-probabilities pbs photon-pola-angle)
      (let* ((p+ (square (cos (- (pbs-angle pbs)
                                 photon-pola-angle))))
             (p- (- 1 p+)))
        (values p+ p-)))

    (define (pbs-tensor pbs photon-pola-angle)
      (let* ((φ₁ (pbs-angle pbs))
             (φ₂ (- π/2 φ₁))
             (θ photon-pola-angle)
             (ampl₁ (abs (cos (- φ₁ θ))))
             (ampl₂ (abs (cos (- φ₂ θ))))
             (θ-string (radians->string θ))
             (φ₁-string (radians->string φ₁))
             (φ₂-string (radians->string φ₂)))
        `((,ampl₁ . ,(string-append θ-string ",+," φ₁-string))
          (,ampl₂ . ,(string-append θ-string ",-," φ₂-string)))))
    
    (define (pbs-activity pbs photon)
      ;; Output (values <photon PBS-ANGLE> #f) if (+) channel.
      ;; Output (values #f <photon PERPENDICULAR-ANGLE>) if (-).
      (define φ (pbs-angle pbs))
      (let-values (((p+ _p-)
                    (pbs-probabilities
                     pbs (photon-polarization-angle photon))))
        (if (< (random-real) p+)
            (values (make-photon φ) #f)
            (values #f (make-photon (+ φ π/2))))))

    (define-record-type <splitter>
      ;; Any of many two-channel devices.  Such a device divides the
      ;; particle beam into proportions cos²(φ) and sin²(φ) if
      ;; behavior is 'optical, or cos²(φ/2) and sin²(φ/2) if behavior
      ;; is 'stern-gerlach. Which particles go which way will be part
      ;; of the design of the device. They are specified by lists and
      ;; distinguished by the ‘equal?’ predicate.
      (%%make-splitter behavior angle cos²-rule sin²-rule)
      splitter?
      (behavior splitter-behavior)
      (angle splitter-angle)
      (cos²-rule splitter-cos²-rule)
      (sin²-rule splitter-sin²-rule))

    (define splitter-cos2-rule splitter-cos²-rule)
    (define splitter-sin2-rule splitter-sin²-rule)

    (define (make-splitter behavior angle cos²-rule sin²-rule)
      (unless (or (eq? behavior 'optical)
                  (eq? behavior 'stern-gerlach))
        (error "make-splitter: expected a 'optical or 'stern-gerlach"
               behavior))
      (unless (real? angle)
        (error "make-splitter: expected a real number" angle))
      (unless (and (pair? cos²-rule)
                   (list? cos²-rule))
        (error "make-splitter: expected a non-empty list" cos²-rule))
      (unless (and (pair? sin²-rule)
                   (list? sin²-rule))
        (error "make-splitter: expected a non-empty list" sin²-rule))
      (%%make-splitter behavior angle cos²-rule sin²-rule))

    (cond-expand
      (chicken (define-record-printer (<splitter> rectype port)
                 (display "<splitter " port)
                 (write (splitter-behavior rectype) port)
                 (display " " port)
                 (write-angle (splitter-angle rectype) port)
                 (display " " port)
                 (write (splitter-cos²-rule rectype) port)
                 (display " " port)
                 (write (splitter-sin²-rule rectype) port)
                 (display ">" port)))
      (else))

    (define (splitter-probabilities splitter particle-kind)
      (let* ((φ (splitter-angle splitter))
             (φ (if (eq? (splitter-behavior splitter) 'stern-gerlach)
                    (* 0.5 φ)
                    φ)))
        (cond ((member particle-kind (splitter-cos²-rule splitter))
               (let* ((p+ (square (cos φ)))
                      (p- (- 1 p+)))
                 (values p+ p-)))
              ((member particle-kind (splitter-sin²-rule splitter))
               (let* ((p+ (square (sin φ)))
                      (p- (- 1 p+)))
                 (values p+ p-)))
              (else (error
                     "splitter-probabilities: not a recognized particle kind"
                     particle-kind)))))

    (define (splitter-activity splitter particle-kind)
      ;; Let the particle pass through into an appropriate channel.
      (let-values (((p+ _p-)
                    (splitter-probabilities splitter particle-kind)))
        (if (< (random-real) p+)
            (values particle-kind #f)
            (values #f particle-kind))))

    (define (estimate-complementary-pair-correlation
             σ_cosφ₁ σ_sinφ₁ σ_cosφ₂ σ_sinφ₂ ; σ = ‘sign’
             fH+V+ fH+V- fH-V+ fH-V-         ; f = ‘frequency’
             fV+H+ fV+H- fV-H+ fV-H-)        ; H = ‘horiz’, V = ‘vert’

      ;;
      ;; Compute an estimate of the correlation coefficient
      ;; -(cos²(φ₁-φ₂)-sin²(φ₁-φ₂)) = -cos(2(φ₁-φ₂) for a
      ;; complementary pair of particles, coherent waves, etc.
      ;;
      ;; This method for estimating correlations is due to
      ;; A. F. Kracklauer.
      ;;

      ;; Compute estimates of products of squares of cosines and sines.
      (let ((cos²φ₁sin²φ₂ (+ fH+V+ fV-H-))
            (cos²φ₁cos²φ₂ (+ fH+V- fV-H+))
            (sin²φ₁sin²φ₂ (+ fH-V+ fV+H-))
            (sin²φ₁cos²φ₂ (+ fH-V- fV+H+)))
        ;; Take square roots. Sometimes one has to use the NEGATIVE
        ;; square root. It depends on the quadrants of φ₁ or φ₂.
        (let ((cosφ₁sinφ₂ (* σ_cosφ₁ σ_sinφ₂ (sqrt cos²φ₁sin²φ₂)))
              (cosφ₁cosφ₂ (* σ_cosφ₁ σ_cosφ₂ (sqrt cos²φ₁cos²φ₂)))
              (sinφ₁sinφ₂ (* σ_sinφ₁ σ_sinφ₂ (sqrt sin²φ₁sin²φ₂)))
              (sinφ₁cosφ₂ (* σ_sinφ₁ σ_cosφ₂ (sqrt sin²φ₁cos²φ₂))))
          ;; Use angle-difference identities. See, for instance, the CRC
          ;; Handbook of Mathematical Sciences, 6th edition, page 170.
          (let ((sin<φ₁-φ₂> (- sinφ₁cosφ₂ cosφ₁sinφ₂))
                (cos<φ₁-φ₂> (+ cosφ₁cosφ₂ sinφ₁sinφ₂)))
            ;; That is it. We have everything we need.
            (- (square sin<φ₁-φ₂>) (square cos<φ₁-φ₂>))))))

    (define (estimate-matched-pair-correlation
             σ_cosφ₁ σ_sinφ₁ σ_cosφ₂ σ_sinφ₂ ; σ = ‘sign’
             fH+V+ fH+V- fH-V+ fH-V-         ; f = ‘frequency’
             fV+H+ fV+H- fV-H+ fV-H-)        ; H = ‘horiz’, V = ‘vert’

      ;; Compute an estimate of the correlation coefficient
      ;; +(cos²(φ₁-φ₂)-sin²(φ₁-φ₂)) = +cos(2(φ₁-φ₂) for a
      ;; matched pair of particles, coherent waves, etc.

      (- (estimate-complementary-pair-correlation
          σ_cosφ₁ σ_sinφ₁ σ_cosφ₂ σ_sinφ₂
          fH+V+ fH+V- fH-V+ fH-V- fV+H+ fV+H- fV-H+ fV-H-)))

    (define (estimate-pair-correlation φ₁ φ₂
                                       optical-or-stern-gerlach
                                       matched-or-complementary
                                       detection-freqs)
      ;;
      ;; Use detection frequencies and trigonometry to get an
      ;; EMPIRICAL ESTIMATE of the value of the correlation.
      ;;
      ;; The ACTUAL values of φ₁ and φ₂ are needed only to determine
      ;; whether a positive or a negative square root is needed in
      ;; certain steps. The signs of the cosine and sine vary
      ;; according to quadrant, but this information is not available
      ;; from the squares of those functions. We have to get the
      ;; information by using the actual angles. However, the actual
      ;; angles ARE available and so it is valid to use them.
      ;;
      (unless (real? φ₁)
        (error "estimate-pair-correlation: expected a real number" φ₁))
      (unless (real? φ₂)
        (error "estimate-pair-correlation: expected a real number" φ₂))
      (unless (or (eq? optical-or-stern-gerlach 'optical)
                  (eq? optical-or-stern-gerlach 'stern-gerlach))
        (error "estimate-pair-correlation: expected 'optical or 'stern-gerlach"
               optical-or-stern-gerlach))
      (unless (or (eq? matched-or-complementary 'matched)
                  (eq? matched-or-complementary 'complementary))
        (error "estimate-pair-correlation: expected 'matched or 'complementary"
               matched-or-complementary))
      (unless (and (list? detection-freqs)
                   (= (length detection-freqs) 8)
                   (every (lambda (x)
                            (and (real? x) (<= 0 x) (<= x 1)))
                          detection-freqs))
        (error "estimate-pair-correlation: expected a list of eight frequencies"
               detection-freqs))

      (let-values (((φ₁ φ₂)
                    (if (eq? optical-or-stern-gerlach 'stern-gerlach)
                        (values (* 1/2 φ₁) (* 1/2 φ₂))
                        (values φ₁ φ₂))))
        (let ((σ_cosφ₁ (if (negative? (cos φ₁)) -1 1))
              (σ_sinφ₁ (if (negative? (sin φ₁)) -1 1))
              (σ_cosφ₂ (if (negative? (cos φ₂)) -1 1))
              (σ_sinφ₂ (if (negative? (sin φ₂)) -1 1)))
          (if (eq? matched-or-complementary 'matched)
              (apply estimate-matched-pair-correlation
                     `(,σ_cosφ₁ ,σ_sinφ₁ ,σ_cosφ₂ ,σ_sinφ₂
                                . ,detection-freqs))
              (apply estimate-complementary-pair-correlation
                     `(,σ_cosφ₁ ,σ_sinφ₁ ,σ_cosφ₂ ,σ_sinφ₂
                                . ,detection-freqs))))))

    (define check-probabilities
      (case-lambda
        ((probs-list) (check-probabilities probs-list 1000))
        ((probs-list
          tolerance) ;; In multiples of one half of fl-epsilon.
         (let ((tol (* tolerance fl-epsilon)))
           (for-each (lambda (p)
                       (unless (and (<= 0 p) (<= (- (+ p p) 2) tol))
                         (error "check-probabilities: probability outside range"
                                p)))
                     probs-list)
           (let ((sum (fold + 0 probs-list)))
             (unless (<= (- (+ sum sum) 2) tol)
               (error "check-probabilities: probabilities sum outside range"
                      probs-list)))))))

    )) ;; end library (epr-simulation)
