# Transients
**Transients in Hippocampal Attractor Networks** - Thesis Project, M.Sc. in Physics of Complex Systems (University of Turin)| *Visiting Research Student* at the University of Ottawa.
Under the supervision of André Longtin (University of Ottawa) and Lamberto Rondoni (University of Turin & Polytechnic University of Turin

## Intro

Computing with a dynamical system implies that this system changes its behavior depending on the quality and quantity of the incoming information.
It is a well-known fact that dissipative networks with symmetric interactions between nodes, operate in a convergent mode, de facto possessing many stable attractors with their basins of attraction likely covering the state space. Such systems can be good models for an associative memory (Hopfield ...). However, computing with attractors makes a very limited use of complex dynamical networks: once an attractor is reached, the "dynamical" nature of the system would no longer be relevant, since wrt the computation, the stytem would now be at a terminal state. The "symmetry assumption" is not only unrealistic in general, but also rule out many possible phenomena and hides many possibilities in their modeling with dynamical systems. 

*Rabinovich et al.* focused they attention on the concept of Winnerless Competition (WLC). WLC itself has been observed in many experiments in hydrodinamics (*Kuppers-Lortz instability* for a large Prandtle number), population biology (*May & Leonard, 1975*) etc.. 

This process can be described by Lotka-Volterra equations (there's a direct connection between Generalized Lotka-Volterra and Rate Models in Neuroscience, see *Rabinovich et al. 2006*), that are well known in population biology and can explain the competition among three species: 

<img src="https://raw.githubusercontent.com/sazio/Transients/master/img/ML3D.png" width="200">



The variables x, y, z are the number of individuals for each population at time *t*, and <img src="https://render.githubusercontent.com/render/math?math=\rho_{ij}"> are competition coefficients measuring how much the *jth* species affects the growth rate of the *ith* species. The non-symmetry of the <img src="https://render.githubusercontent.com/render/math?math=\rho_{ij}"> guarantees the WLC behavior of the above mentioned dynamical system. The phase portrait of such a behavior is a heteroclinic contour (see below), there's a sequential switching from saddle point to saddle point.


<img src="https://raw.githubusercontent.com/sazio/Transients/master/img/phase_portraitml.png" width="600">

And what about the solutions? You can clearly see a sequential switching 

<img src="https://raw.githubusercontent.com/sazio/Transients/master/img/solutions.png" width="600">

In *Afraimovich et al. 2003* conditions for the robustness of the WLC, existence and stability of the heteroclinic contour are further investigated. 

I've also tried to reproduce a WLC competition with a network of 3 FitzHugh-Nagumo neurons (when N > 3, the dynamics of the system can be very complex and even chaotic *Varona et al. 2001* --> WLC in molluscan hunting behaviour) as in *Rabinovich et al. 2001*. 

(Unusual FHN model, with 3 equations instead of 2) 


<img src="https://raw.githubusercontent.com/sazio/Transients/master/img/fhn3eqs.png" width="400">


This is the projection of a 9-dimensional heteroclinic orbit of three inhibitory coupled FHN neirons in a 3 Dim space, <img src="https://render.githubusercontent.com/render/math?math=\xi_{i}">s are linear combinations of the actual phase space variables of the system. 
<img src="https://raw.githubusercontent.com/sazio/Transients/master/img/FHN_9D.png" width="600">


I think I've gotten something similar - not that cute but ok-, still trying to understand how Julia's plot works but we have 3 "limit" cycles... (there are a few theorems on the birth of a stable limit cycle in the case on an appropriate perturbation, I'm looking through all of those right now)

<img src="https://raw.githubusercontent.com/sazio/Transients/master/img/phase_portraitFHN.png" width="600">


## Some words on Winnerless Competition

Competition without a winner is a widely known phenomenon in systems involving more than two interacting agents that satisfies a relationship similar to the voting paradox (Game Theory involved) or the popular game, *rock-paper-scissors*. Such interactions lead to a *nontransitive competition* i.e. a cyclic behavior. The mathematical image of this is a robust heteroclinic cycle which was first formulated in *Busse & Heikes, 1980* in modeling convection in a rotating layer. 

<img src="https://raw.githubusercontent.com/sazio/Transients/master/img/SHC.png" width="600">

A segment of a stable heteroclinic chain adopted from [Afraimovich et al., 2004]. Each node in the chain represents a saddle possessing both stable and unstable manifolds. Note that the trajectories do not necessarily lead to the saddle node that is being approached; some may slowly pass through its vicinity as the blue and the green paths do. The distance from the state vector to the saddle is proportional to the speed. The chain itself is stable, so that the blue and the green paths wandering in the saddle neighborhoods are not limit cycles; they can be any solution staying in the vicinity of the chain.


<img src="" width="600">


<img src="" width="600">
<img src="" width="600">

