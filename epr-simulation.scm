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

  (export π/180 π/8 π/4 π3/8 π/2 π)

  (export degrees->radians radians->degrees)

  (export <photon>
          ;; A plane-polarized photon.
          make-photon
          photon?
          photon-polarization-angle
          *photon-pair-probability* ; Parameter, default 0.5
          photon-pair-probabilities ; Two values, default 0.5, 0.5
          photon-pair-source)       ; Two photons, complementary pair.

  (export <pbs>
          ;; Polarizing beam-splitter with plane-polarized output.
          make-pbs
          pbs?
          pbs-angle-in      ; Angle wrt the incident photon.
          pbs-angle-out+    ; Output angle wrt the next arm, or #f.
          pbs-angle-out-    ; Output angle wrt to the next arm, or #f.
          pbs-probabilities ; Probabilities given an incident photon.
          pbs-activity)     ; Activity given an incident photon.

    (export <sgm>
          ;; A two-channel device that behaves like a Stern-Gerlach
          ;; magnet.
          make-sgm
          sgm?
          sgm-angle-in                ; Angle wrt the incident object.
          sgm-cos²-rule-particles     ; Particles that obey cos² rule.
          sgm-cos2-rule-particles     ; Synonym.
          sgm-sin²-rule-particles     ; Particles that obey sin² rule.
          sgm-sin2-rule-particles     ; Synonym.
          sgm-probabilities  ; Probabilities given an incident object.
          sgm-activity)      ; Activity given an incident object.

  (import (scheme base)
          (scheme case-lambda)
          (scheme inexact)
          (only (srfi 27) random-real)
          (only (srfi 144) fl-epsilon))

  (cond-expand
    (chicken (import (only (chicken base) define-record-printer)
                     (scheme write)))
    (else))

  (begin

    (cond-expand
      (chicken (define (write-angle angle port)
                 (let* ((angle/π (/ angle π))
                        (angle/π*64 (* 64 angle/π))
                        (iangle/π*64 (round angle/π*64))
                        (diff (abs (- angle/π*64 iangle/π*64)))
                        (exact-enough (<= diff (* 500 fl-epsilon
                                                  (abs angle/π*64)))))
                   (display "π×" port)
                   (let ((angle/π
                          (if exact-enough
                              (exact angle/π)
                              angle/π)))
                     (write angle/π port)))))
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

    ;; OEIS A000796
    (define π 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214)

    (define (degrees->radians x) (* x π/180))
    (define (radians->degrees x) (/ x π/180))

    (define-record-type <photon>
      ;; A photon is a light pulse of unit intensity. It has some
      ;; polarization angle orthogonal to its direction of motion. The
      ;; polarization angle is given by a real number.
      (%%make-photon angle)
      photon?
      (angle photon-polarization-angle))

    (define (make-photon angle)
      (unless (real? angle)
        (error "make-photon: expected a real number" angle))
      (%%make-photon angle))

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
      (%%make-pbs angle-in angle-out+ angle-out-)
      pbs?
      (angle-in pbs-angle-in)
      (angle-out+ pbs-angle-out+)
      (angle-out- pbs-angle-out-))

    (define make-pbs
      (case-lambda
        ((angle-in angle-out+ angle-out-)
         (unless (real? angle-in)
           (error "make-pbs: expected a real number" angle-in))
         (unless (or (not angle-out+) (real? angle-out+))
           (error "make-pbs: expected a real number or #f"
                  angle-out+))
         (unless (or (not angle-out-) (real? angle-out-))
           (error "make-pbs: expected a real number or #f"
                  angle-out-))
         (%%make-pbs angle-in angle-out+ angle-out-))
        ((angle-in)
         (make-pbs angle-in #f #f))))

    (cond-expand
      (chicken (define-record-printer (<pbs> rectype port)
                 (display "<pbs " port)
                 (write-angle (pbs-angle-in rectype) port)
                 (display " " port)
                 (let ((x (pbs-angle-out+ rectype)))
                   (if (real? x)
                       (write-angle x port)
                       (write x port)))
                 (display " " port)
                 (let ((x (pbs-angle-out- rectype)))
                   (if (real? x)
                       (write-angle x port)
                       (write x port)))
                 (display ">" port)))
      (else))

    (define (pbs-probabilities pbs photon-pola-angle)
      (let* ((p+ (square (cos (- (pbs-angle-in pbs)
                                 photon-pola-angle))))
             (p- (- 1 p+)))
        (values p+ p-)))

    (define (pbs-activity pbs photon)
      ;; Output (values #t #f) or (values <photon ANGLE-OUT> #f)
      ;;   if (+) channel is chosen.
      ;; Output (values #f #t) or (values #f <photon ANGLE-OUT>)
      ;;   if (-) channel is chosen.
      (let-values (((p+ _p-)
                    (pbs-probabilities
                     pbs (photon-polarization-angle photon))))
        (if (< (random-real) p+)
            (let ((angle-out (pbs-angle-out+ pbs)))
              (values (if angle-out (make-photon angle-out) #t)
                      #f))
            (let ((angle-out (pbs-angle-out- pbs)))
              (values #f
                      (if angle-out (make-photon angle-out) #t))))))

    (define-record-type <sgm>
      ;; Any two-channel device similar to a Stern-Gerlach magnet.
      ;; Such a device divides the particle beam into proportions
      ;; cos²(φ/2) and sin²(φ/2). Which particles go which way will be
      ;; part of the design of the device. They are specified by lists
      ;; and distinguished by the ‘equal?’ predicate.
      (%%make-sgm angle-in cos²-rule-particles sin²-rule-particles)
      sgm?
      (angle-in sgm-angle-in)
      (cos²-rule-particles sgm-cos²-rule-particles)
      (sin²-rule-particles sgm-sin²-rule-particles))

    (define sgm-cos2-rule-particles sgm-cos²-rule-particles)
    (define sgm-sin2-rule-particles sgm-sin²-rule-particles)

    (define (make-sgm angle-in
                      cos²-rule-particles sin²-rule-particles)
      (unless (real? angle-in)
        (error "make-sgm: expected a real number" angle-in))
      (unless (and (pair? cos²-rule-particles)
                   (list? cos²-rule-particles))
        (error "make-sgm: expected a non-empty list"
               cos²-rule-particles))
      (unless (and (pair? sin²-rule-particles)
                   (list? sin²-rule-particles))
        (error "make-sgm: expected a non-empty list"
               sin²-rule-particles))
      (%%make-sgm angle-in cos²-rule-particles sin²-rule-particles))

    (cond-expand
      (chicken (define-record-printer (<sgm> rectype port)
                 (display "<sgm " port)
                 (write-angle (sgm-angle-in rectype) port)
                 (display " " port)
                 (write (sgm-cos²-rule-particles rectype) port)
                 (display " " port)
                 (write (sgm-sin²-rule-particles rectype) port)
                 (display ">" port)))
      (else))

    (define (sgm-probabilities sgm particle-kind)
      (define φ/2 (* 0.5 (sgm-angle-in sgm)))
      (cond ((member particle-kind (sgm-cos²-rule-particles sgm))
             (let* ((p+ (square (cos φ/2)))
                    (p- (- 1 p+)))
               (values p+ p-)))
            ((member particle-kind (sgm-sin²-rule-particles sgm))
             (let* ((p+ (square (sin φ/2)))
                    (p- (- 1 p+)))
               (values p+ p-)))
            (else (error
                   "sgm-probabilities: not a recognized particle kind"
                   particle-kind))))

    (define (sgm-activity sgm particle-kind)
      ;; Let the particle pass through into an appropriate channel.
      (let-values (((p+ _p-) (sgm-probabilities sgm particle-kind)))
        (if (< (random-real) p+)
            (values particle-kind #f)
            (values #f particle-kind))))

    )) ;; end library (epr-simulation)
