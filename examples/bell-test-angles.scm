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

(import (scheme base)
        (scheme write)
        (epr-simulation))

(cond-expand
  (chicken (import (format)))
  (gauche (import (srfi 28))))    ; FIXME: not tested yet with Gauche.

(define bell-test-angles
  `((0 π/8)
    (0 π3/8)
    (π/4 π/8)
    (π/4 π3/8)))

(define (simulate-one-event φᴸ φᴿ)
  (let ((pbsᴸ (make-pbs φᴸ))
        (pbsᴿ (make-pbs φᴿ)))
    (let*-values (((Pᴸ Pᴿ) (photon-pair-probabilities))
                  ((Aᴸ Aᴿ) (photon-pair-amplitudes))
                  ((ξᴸ ξᴿ) (photon-pair-source 0 π/2))
                  ((detectedᴸ+ detectedᴸ-) (pbs-activity pbsᴸ ξᴸ))
                  ((detectedᴿ+ detectedᴿ-) (pbs-activity pbsᴿ ξᴿ)))
      `(,(if detectedᴸ+ '+ '-) ,(if detectedᴿ+ '+ '-)))))

(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
(write (simulate-one-event 0 π/8))(newline)
