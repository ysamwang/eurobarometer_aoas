The original Stata data is in ZA1961_v1-0-1.dta

The euro34.R file cleans the data. Specifically, it removes individuals who didn't answer all 3 questions or anyone who gave an answer with a tie. This leaves 11872 out of the original 12733 individuals. This file also takes from the raw stata file into the cleanEuro34.txt file.

In the cleanEuro34.txt file, you have the cleaned responses for each individual. Each row corresponds to an individual In each column, you have the ranking given by the individual for a specific policy. So if the 1st policy was the 2nd most prefered policy for individual 100, then the (100, 1) element will be 2. If there is an NA, that means that the individual did not give a rank to that policy. The specific policies are given in the eurobarometerSummary.pdf file (section 1.2)
