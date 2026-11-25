# Fuzzy
Enhancing Particle Swarm Optimization via a Fuzzy Logic Add-On: A Comparative Study across PSO Variants

## Abstract

In Particle Swarm Optimization (PSO), the inertia weight (ω) determines the extent to which a particle's previous velocity is carried over to the subsequent iteration. At the same time, the cognitive coefficient (c1) scales the cognitive component of the velocity update. Most PSO variants assign these parameters fixed values or vary them according to a predefined schedule. This study presents a Mamdani fuzzy controller that adjusts both parameters in response to the current state of the swarm. At each iteration, the controller utilizes a set of if–then rules to determine ω and c1 from two inputs: a stagnation counter and a measure of swarm diversity computed across individual dimensions. We integrated the controller into BLPSO, CLPSO, and PSOGSA to evaluate its use across three distinct learning mechanisms. For each algorithm, we compared the fuzzy variant with its original implementation on 40 benchmark functions in both 30 and 50 dimensions. We assessed differences between each pair using the Wilcoxon rank-sum test on results from 30 runs per configuration.

<img width="3200" height="990" alt="fuzzyflow" src="https://github.com/user-attachments/assets/1830dd3c-899c-44c2-bda8-7f78d0ab34ea" />


