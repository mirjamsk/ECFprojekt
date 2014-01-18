CLONAL algorithm for ECFramework
=============================
/**
 * \brief Clonal Selection Algorithm (see e.g. http://en.wikipedia.org/wiki/Clonal_Selection_Algorithm)
 * 
 * this CLONALG implements: 
 *			- cloning Versions:	- static cloning :  n of the best antibodies are cloned beta*populationSize times
 *						- proportional cloning: number of clones per antibody is proportional to that ab's fitness
 *			- inversely proportional hypermutation : better antibodies are mutated less
 *			- selectionSchemes:	- CLONALG1 - at new generation each antibody will be substituded by the best individual of its set of beta*population clones
 *						- CLONALG2 - new generation will be formed by the best (1-d)*populationSize clones ( or all if the number of clones is less than that )
 *          		- birthPhase: where d * populationSize of new antibodies are randomly created and added to the population for diversification
 *							
 * CLONALG algorithm accepts only a single FloatingPoint genotype
 * Additionally, if chosen, selectionScheme CLONALG2 adds a FloatingPoint genotype  (parentAntibody) to mark which clone came from which antibods 
 */
=============================
important: main.cpp is from ECF_1.3/examples/COCO/
           it does not implement AGING! 
