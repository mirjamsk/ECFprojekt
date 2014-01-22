ABC algorithm for ECFramework
=============================
/**
 * \brief Artificial Bee Colony algorithm (see e.g.http://www.scholarpedia.org/article/Artificial_bee_colony_algorithm)
 * 
 * ABC algorithm accepts only a single FloatingPoint genotype (vector of real values).
 * This version of ABC algorithm uses a built-in ECF operator for selecting the best individuals
 * Additionally, it adds the following genotype for algorithm implementation:
 * 		- trial: floatingPoint genotype serving as a generation counter for each individual
 */


=============================
important: main.cpp is from ECF_1.3/examples/COCO/