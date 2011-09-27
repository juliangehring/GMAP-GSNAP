#ifndef MAXENT_INCLUDED
#define MAXENT_INCLUDED

#define DONOR_MODEL_LEFT_MARGIN 3 /* Amount in exon.  Does not include GT */
#define DONOR_MODEL_RIGHT_MARGIN 6 /* Amount in intron */

#define ACCEPTOR_MODEL_LEFT_MARGIN 20 /* Amount in intron.  Includes AG */
#define ACCEPTOR_MODEL_RIGHT_MARGIN 3 /* Amount in exon */

extern double
Maxent_donor_prob (char *sequence);

extern double
Maxent_acceptor_prob (char *sequence);

#endif

