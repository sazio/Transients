# Transients
**Transients in Hippocampal Attractor Networks** - Thesis Project, M.Sc. in Physics of Complex Systems (University of Turin)| *Visiting Research Student* at the University of Ottawa.
Under the supervision of André Longtin (University of Ottawa) and Lamberto Rondoni (University of Turin & Polytechnic University of Turin)

## Intro

Computing with a dynamical system implies that this system changes its behavior depending on the quality and quantity of the incoming information.
It is a well-known fact that dissipative networks with symmetric interactions between nodes, operate in a convergent mode, de facto possessing many stable attractors with their basins of attraction likely covering the state space. Such systems can be good models for an associative memory (Hopfield ...). However, computing with attractors makes a very limited use of complex dynamical networks: once an attractor is reached, the "dynamical" nature of the system would no longer be relevant, since wrt the computation, the stytem would now be at a terminal state. The "symmetry assumption" is not only unrealistic in general, but also rule out many possible phenomena and hides many possibilities in their modeling with dynamical systems. 

*Rabinovich et al.* focused their attention on the concept of Winnerless Competition (WLC). WLC itself has been observed in many experiments in hydrodinamics (*Kuppers-Lortz instability* for a large Prandtle number), population biology (*May & Leonard, 1975*) etc.. 

This process can be described by Lotka-Volterra equations (there's a direct connection between Generalized Lotka-Volterra and Rate Models in Neuroscience, see *Rabinovich et al. 2006*, *Fukai & Tanaka, 1997*), that are well known in population biology and can explain the competition among three species: 

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


## Some words on the Winnerless Competition

Competition without a winner is a widely known phenomenon in systems involving more than two interacting agents that satisfies a relationship similar to the voting paradox (Game Theory involved) or the popular game, *rock-paper-scissors*. Such interactions lead to a *nontransitive competition* i.e. a cyclic behavior. The mathematical image of this is a robust heteroclinic cycle which was first formulated in *Busse & Heikes, 1980* in modeling convection in a rotating layer. 

<img src="https://raw.githubusercontent.com/sazio/Transients/master/img/SHC.png" width="600">

A segment of a stable heteroclinic chain adopted from [Afraimovich et al., 2004]. Each node in the chain represents a saddle possessing both stable and unstable manifolds. Note that the trajectories do not necessarily lead to the saddle node that is being approached; some may slowly pass through its vicinity as the blue and the green paths do. The distance from the state vector to the saddle is proportional to the speed. The chain itself is stable, so that the blue and the green paths wandering in the saddle neighborhoods are not limit cycles; they can be any solution staying in the vicinity of the chain.


## Unconscious Animal Behaviour

Reproducible sequences are universal objects in unconscious animal behavior. In the course of evolution, sequences are optimized in some sense and became generic for different animals. The search for the food, for example, is a key task for survival and different animals demonstrate pro- totype patterns in the search for their preys. These patterns look like a prescribed combination of rotation and long-range trajectories in a sequential manner. *Dynamical activity of net- works consisting of competitive neurons can explain the emergence of such patterns.* 

Varona et al. 2002, have been working on *Clione Limacina*, a blind molluscan, which possesses a gravity sensory organ, called stratocyst, for orientation. The stratocyst contains a free stone-like mass, called stratolyth, which signals the direction of motion by exciting sensory neurons around it.

At a given time, the animal's locomotion is in one of the two characteristic models: navigation and hunting. 
In navigation mode, clione's CPG (Central Pattern generator) produces a uniform signal for the actuators (wings and tail) depending on the target's location and the actual direction. In the absence of strong perturbations, such a contro system should generate a uniform control signal (see below). 

On the other hand there's the hunting mode, which triggers a complex dynamics in the pattern generation which produces a sequential control signal for the actuators in order to approach the prey in a non-uniform manner. This sequence, however, cannot be a totally random (or chaotic) signal as such control signal would cause some sort of a *Brownian motion*. In fact, the animal has been shown to follow a prescribed strategy, i.e. a reproducible sequence of switchings, in the hunting mode which can be explained by the **winnerless competition** framework.


<img src="https://github.com/sazio/Transients/blob/master/img/Stratolyth.png" width="600">


<img src="https://github.com/sazio/Transients/blob/master/img/sequence.png" width="600">

In the navigation mode, the statolyth excites one sensory neuron (the unique winner, winner-take-all), depending on the moving direction. Due to the inhibitory connections in the sensory network of the statocyst, all other neurons are kept silent (see the first sequence). The winner neuron signals the actual direction of the animal. The brain evaluates this feedback and actuates the wings and tail to maintain the desired movement direction. On the right, in the hunting mode (triggered by an auxiliary signal from the H (hunting) neuron), the inhibitory interaction among sensory neurons conducts a winnerless competition, where sensory neurons fire in a sequence. The switching times may be irregular, but the sequence of the winners is deterministic and reproducible (see the second sequence). When hunting, the effect of the statolyth is ignored, and the robust sequence, as prescribed by the nature in evolution, results in a sequence of actions on the tail and the wings.

## Binding and Chunking Dynamics (Kiebel and Friston 2012, Rabinovich et al 2014)

What if we'd like to generate "non-repetitive" sequences? A Stable Heteroclinic Channel (SHC) can only generate stereotyped sequences. For example, in the previous case where we had 3 saddle points (GLV), a SHC forces the trajectories trough the saddle points in a very regular manner, e.g. 1-2-3-1-2-3.... In contrast, it can't generate something like "1-2-3-2-1-3 ... " because in this case the sequence is not repetitive. However, to model sensory input, we'd like to be able to recombine basic elements of this sequence in ever-changing sequences. This concept can be easily illustrated on a language example (*Cona and Semenza, 2017*), since, sequence of letters compound syllables, syllables compound words, words compound sentences and so on. Language is , in fact, a hierarchical sequential process. 

A plausible solution might be to construct a hierarchy of SHCs, which can encode sequences generated by SHCs whose attractor topology is changed by a supraordinate SHC. This can be achieved by making the connectivity matrix at a subordinate level a function of the output states of the supraordinate level. This enables the hierarchy to generate sequences of sequences to any hierarchical depth required (i.e. a Multilevel Network). 

This hierarchical modelling framework can be helpful for understanding chunking phenomena. Chunking is a dynamical phenomenon that nature uses to perform information processing of long sequences by dividing them in shorter information items. It is pretty useful for example, in making an efficient use of short-term memory by breaking up long strings of information (e.g. in language, one can see the separation of a novel on chapters, paragraphs, sentences and words...) 

In Rabinovich et al 2014, a cognitive network architecture that hierarchically chunks and super-chunks switching sequences of metastable states produce by WLC heteroclinic dynamics has been proposed. 

The concept of chunk has been firstly introduced by *Miller 1956*. The key notion here, was that short-term storaging is not rigid at all but open to strategies that can expand its capacities (e.g. chunking). 

Chunking involves two processes: 1) concatenation of units in a block and 2) segmentation of the blocks. As we've previously said, it is related to the hierarchical organization of perceptual, cognitive or behavioural sequential activity. In motor control (*Rosenbaum et al, 1983*) sequences can consist of subsequences and these can in turn consist of sub-sub-sequences. *Clerget et al. 2013* hypothesize that motor behavior shares some similarities with language, namely that a complex action can be viewed as a chain of subordinate movements, which need to be combined according to certain rules in order to reach a given goal. 

Here metastability is a key element of transient cognitive dynamics participating in chunking processes. In particular, such dynamics can be represented as a sequential switching between different metastable states (*Rabinovich et al, 2008a,b*).  


## References 
Clerget, E., Andres, M., and Olivier, E. (2013). Deficit in complex sequence processing after a virtual lesion of left BA45. PLoS ONE 8:e63722. doi: 10.1371/journal.pone.0063722

Kiebel, S. J., Friston, K. J. (2012), Recognition of Sequences of Sequences Using Nonlinear Dynamical Systems, Chapter 6 of Principles of Brain Dynamics (Rabinovich, Varona), doi: 10.7551/mitpress/9108.001.0001.

Miller, G. A. (1956). The magical number seven plus or minus two: some limits on our capacity for processing information. Psychol. Rev. 63, 81–97. doi: 10.1037/h0043158

Rosenbaum, D. A., Kenny, S. B., and Derr, M. A. (1983). Hierarchical control of rapid movement sequences. J. Exp. Psychol. Hum. Percept. Perform. 9, 86–102. doi: 10.1037/0096-1523.9.1.86

Rabinovich, M. I., Varona, P., Tristan, I., Afraimovich, V. S. (2014), Chunking dynamics: heteroclinics in mind, Front. Comput. Neurosci.,  doi: 10.3389/fncom.2014.00022

Rabinovich, M., Huerta, R., and Laurent, G. (2008a). Neuroscience. Transient dynamics for neural processing. Science 321, 48–50. doi: 10.1126/science. 1155564

Rabinovich, M. I., Huerta, R., Varona, P., and Afraimovich, V. S. (2008b). Transient cognitive dynamics, metastability, and decision making. PLoS Comput Biol 4:e1000072. doi: 10.1371/journal.pcbi.1000072




