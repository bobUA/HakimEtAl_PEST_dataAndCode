Data and code to generate Figures from 

Hakim et al. Evaluating the cognitive mechanisms of phishing detection with PEST, an ecologically valid lab-based measure of phishing susceptibility

NOTE: Figure 4 requires data from the original PHIT task. These are available 

Data files are csv files.  Naming has the following form:

	scamdata_SUBJECTNUMBER_DATETIME_AGE_GENDER.dat

		e.g. scamdata_1_10Oct2018090103_18_F.dat

Each datafile has 7 columns :

userId    : subject response (1 - safe with high confidence, 2 - safe with low confidence, 3 - scam with low confidence, 4 - scam with high confidence)
reactTime : reaction time in seconds 
category  : PHIT Email Category (and custom categories for pooled scam/safe emails)
type      : weapon of influence (for PHIT emails only)
hasAtt    : binary indicating whether email has an attachment
realID    : real email identifier (scam or safe)
emailCode : unique ID of each email - used to locate specific emails within excel files

Place to start with analysis is 

	main.m

This will recreate all figures and analyses from the paper (except the task illustration, Figure 1)