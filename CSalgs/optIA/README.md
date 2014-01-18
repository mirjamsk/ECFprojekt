Optimization Immune  Algorithm (opt-IA) for ECFramework
---


+ This opt-IA implements:  
 + static cloning: all antibodies are cloned _dup_ times, making the size of the clone population equal _dup*spoplationSize_
 + inversely proportional hypermutation: better antibodies are mutated less
 + static pure aging: if an antibody exceeds tauB number of trials, it is replaced with a new randomly created antibody
 + birthPhase: if the number of antibodies that survive the aging Phase is less than populationSize, new randomly created abs are added to the population
 + optional elitism


+ opt-IA algorithm accepts only a single FloatingPoint genotype

+ Additionally, opt-IA adds a FloatingPoint genotype (age) 
 
=============================


*important: main.cpp is from ECF_1.3/examples/COCO/*
