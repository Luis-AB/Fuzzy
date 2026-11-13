# Fuzzy
Enhancing Particle Swarm Optimization via a Fuzzy Logic Add-On: A Comparative Study across PSO Variants

## Abstract

Particle Swarm Optimization (PSO) is a popular metaheuristic inspired by collective behavior, where individuals learn from their own experience and the success of others. Its performance depends on the inertia weight (ω) and the acceleration coefficients, which most variants fix or schedule in advance.

This work proposes a Mamdani fuzzy controller, which expresses the adjustment as if-then rules, as an add-on that sets ω and the cognitive coefficient (c₁) at each iteration based on two indicators of the swarm state: a stagnation counter and dimension-wise diversity. We attached it to three PSO variants with different learning mechanisms (BLPSO, CLPSO, and PSOGSA) and compared each with its original version on 40 benchmark functions in 30 and 50 dimensions, using the Wilcoxon rank-sum test over 30 runs.


