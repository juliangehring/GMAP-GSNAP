/* $Id: complement.h,v 1.6 2005/02/15 01:50:55 twu Exp $ */
#ifndef COMPLEMENT_INCLUDED
#define COMPLEMENT_INCLUDED

/* Purine (A,G): R => T,C: Y
   Pyrmidine (C,T): Y => G,A: R
   A,T: W => T,A: W
   G,C: S => C,G: S
   A,C: M => T,G: K
   G,T: K => C,A: M
   A,T,C: H => T,A,G: D
   G,C,T: B => C,G,A: V
   G,A,C: V => C,T,G: B
   G,A,T: D => C,T,A: H
   G,A,T,C: N => C,T,A,G: N 

   Also X
*/

/* Also provides reverse of comp symbols |, >, <, ), (, ], [, =, # */
   
/*
                              1111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990000000000111111111122222222
                    01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567
                    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~@
*/

#define COMPLEMENT "???????????????????????????????? ??#$%&')(*+,-./0123456789:;>=<??TVGHEFCDIJMLKNOPQYSAUBWXRZ]?[^_`tvghefcdijmlknopqysaubwxrz}|{~?"

#endif
