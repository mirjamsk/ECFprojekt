CLONAL algorithm for ECFramework
=============================
 * \brief Clonal Selection Algorithm (see e.g. http://en.wikipedia.org/wiki/Clonal_Selection_Algorithm)
 * 
 * this CLONALG implements:
	- static cloning :  n of the best antibodies are cloned beta time, making the size of the clones population  equal n*beta
 	- inversely proportional hypermutation : better antibodies are mutated less
 	- CLONALG2 - keeps best (1-d)*populationSize antibodies ( or all if the number of clones is less than that )
	- birthPhase where d * populationSize of new antibodies are randomly created and added to the population
 	
 *
 * CLONALG algorithm accepts only a single FloatingPoint genotype (vector of real values).
 * 
 */

=============================
important: main.cpp is from ECF_1.3/examples/COCO/
it does not implement AGING!