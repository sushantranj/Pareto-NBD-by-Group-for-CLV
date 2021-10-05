# Pareto-NBD-by-Group-for-CLV
Implementation of Pareto NBD by group in R

The Pareto/NBD model makes the following assumptions regarding the customer population: 
•	Purchase count follows a Poisson distribution with rate λ. In other words, the timing of these purchases is somewhat random, but the rate (in counts/unit time) is constant. In turn, this implies that the inter-purchase time at the customer level should follow an exponential distribution.
•	Lifetime distribution follows an exponential distribution with slope μ. The expectation value of such distribution is 1/μ and corresponds to the lifetime of the user. 
•	The latent parameters λ and μ are constrained by two prior gamma distributions representing our belief of how these latent parameters are distributed among the population of customers. These two gamma distributions have parameters (r, α) for the purchase count and (s, β) for the lifetime. The goal is to find these four parameters. From these, all actionable metrics can be derived.


The input data is in the marketing.csv where the variable region is used as a group variable

R supports packages like BTYD and BTYDPlus to implement the probability distribution models related to the customer lifetime value modeling. Because of the manner (full transpose) in which the calibration matrix is generated, these fail for a larger set of data. In order to get the problem solved through R, the core components of matrix creation and maximum likelihood estimation was rewritten as a set of functions ( as mentioned in the research papers mentioned in the appendix)would be in 

The output would be in the below format 

cust	region	  expected_transactions	probability_of_activity 	residual_transactions	expected_purchase_from_a_new_member_bygroup 
a100000004	House Ads	1	0.00	1.00	3
![image](https://user-images.githubusercontent.com/88220240/128811495-60aa3dee-c83b-4a2a-8dba-70dde9196ea9.png)



The above output could be married with the average revenue of the member along side an appropriate discount to get to the final customer lifetime value .
