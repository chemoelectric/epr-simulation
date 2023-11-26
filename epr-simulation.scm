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
          make-photon
          photon?
          photon-polarization-angle)
  (export <pbs>
          ;; Polarizing beam-splitter with plane-polarized output.
          make-pbs
          pbs?
          pbs-angle-in
          pbs-angle-out+
          pbs-angle-out-)

  (import (scheme base)
          (scheme inexact)
          (only (srfi 144) fl-epsilon))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;          (srfi 27))                    ; Random bits.

  (cond-expand
    (chicken (import (only (chicken base) define-record-printer)
                     (scheme write)))
    (else))

  (begin

    (cond-expand
      (chicken (define (write-angle angle port)
                 (let* ((angle/π (/ angle π))
                        (angle/π*8 (* 8 angle/π))
                        (iangle/π*8 (round angle/π*8))
                        (diff (abs (- angle/π*8 iangle/π*8)))
                        (exact-enough (<= diff (* 500 fl-epsilon
                                                  (abs angle/π*8)))))
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

    (define-record-type <pbs>
          ;; Polarizing beam-splitter with plane-polarized output.
      (%%make-pbs angle-in angle-out+ angle-out-)
      pbs?
      (angle-in pbs-angle-in)
      (angle-out+ pbs-angle-out+)
      (angle-out- pbs-angle-out-))

    (define (make-pbs angle-in angle-out+ angle-out-)
      (unless (real? angle-in)
        (error "make-pbs: expected a real number" angle-in))
      (unless (real? angle-out+)
        (error "make-pbs: expected a real number" angle-out+))
      (unless (real? angle-out-)
        (error "make-pbs: expected a real number" angle-out-))
      (%%make-pbs angle-in angle-out+ angle-out-))

    (cond-expand
      (chicken (define-record-printer (<pbs> rectype port)
                 (display "<pbs " port)
                 (write-angle (pbs-angle-in rectype) port)
                 (display " " port)
                 (write-angle (pbs-angle-out+ rectype) port)
                 (display " " port)
                 (write-angle (pbs-angle-out- rectype) port)
                 (display ">" port)))
      (else))

    )) ;; end library (epr-simulation)
