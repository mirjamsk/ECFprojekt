
CLONAL algorithm for ECFramework
===


*(see e.g. http://en.wikipedia.org/wiki/Clonal_Selection_Algorithm/)*


+ This CLONALG implements :
	+ cloning Versions:
	 + static cloning:  n of the best antibodies are cloned _beta*populationSize_ times
	 + proportional cloning:  number of clones per antibody is proportional to that ab's fitness
	 + inversely proportional hypermutation: better antibodies are mutated less
	+ selectionSchemes:
	 + CLONALG1: at new generation each antibody will be substituded by the best individual of its set of _beta*population_ clones
	 + CLONALG2: new generation will be formed by the best _(1-d)*populationSize_ clones ( or all if the number of clones is less than that )
	+ birthPhase: where _d*populationSize_ of new antibodies are randomly created and added to the population for diversification
           
+ CLONALG algorithm accepts only a single FloatingPoint genotype
+ Additionally, if chosen, selectionScheme CLONALG1 adds a FloatingPoint genotype  (parentAntibody) to mark which clone came from which antibody


===


*important: main.cpp is from ECF_1.3/examples/COCO/*
*it does not implement AGING!*
