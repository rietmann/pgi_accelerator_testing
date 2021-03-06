/**
 * Patus Runtime Library
 */

#ifndef __PATUSRT_H__
#define __PATUSRT_H__

#ifdef __cplusplus
extern "C" {
#endif


/*************************************************************************/
/* Timer Functions                                                       */
/*************************************************************************/

/**
 * Starts the timer.
 */
void tic ();

/**
 * Stops the timer and prints a time value to stdout.
 */
void toc (long nFlopsPerStencil, long nStencilComputationsCount, long nBytesTransferred);


#ifdef __cplusplus
}
#endif

#endif
