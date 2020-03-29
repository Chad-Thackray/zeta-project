# zeta-project

This is a collection of tools for verifying that the zeros of the Riemann Zeta function lie on the critical line up to a given gram point. Originally created for a dissertation detailing the mathematics behind such algorithms. 

Zero_confirm_accurate Uses up to the third remainder term and is thus markedly more accurate, which is required when dealing with higher values of t. 

Zero_confirm_first only uses the first remainder term and is therefore significantly faster in verifying the zeros. 


Instructions for use:

Scroll down to the "activation code" section and enter the range of gram points between which you would like the Riemann Hypothesis to be verified. Simply change the range of numbers in the verifyBulk function called immediately after the opening of the if statement. Recall that this program is technically only valid at verification past gram point 300, although it will give the correct answer for all gram points. 

By default we examine the range [g_300, g_99,999], the same as in the original dissertation.

If the program takes an unusually large amount of time when processing the Lehman groups, this usually means that it has encountered an extreme case of Lehmer's phenomenon and requires more accuracy to discern the appropriate value of h_n. 


Tips for optimization:

Depending on the range you've selected and the amount of threads your machine can run, increasing or decreasing the groupsize variable may help.

In the "Error Terms" section one may change the cut_off variable. Increasing the value includes higher degree terms in the series expansion. Increasing accuracy but also taking longer to compute


I have only tested these algorithms up to gram point 250,000. I cannot vouch for their behaivior beyond that. It is likely some functions would need to be overhauled and perhaps faster computational methods considered. 



