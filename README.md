# Transients
**Transients in Hippocampal Attractor Networks** - Thesis Project, M.Sc. in Physics of Complex Systems (University of Turin)| *Visiting Research Student* at the University of Ottawa.
Under the supervision of André Longtin (University of Ottawa) and Lamberto Rondoni (University of Turin & Polytechnic University of Turin

## Intro

Computing with a dynamical system implies that this system changes its behavior depending on the quality and quantity of the incoming information. It is an enormous field, *Rabinovich et al.* focused they attention on the concept of Winnerless Competition (WLC). WLC itself has been observed in many experiments in hydrodinamics (*Kuppers-Lortz instability* for a large Prandtle number), population biology (*May & Leonard, 1975*) etc.. 

This process can be described by Lotka-Volterra equations (there's a direct connection between Generalized Lotka-Volterra and Rate Models in Neuroscience, see *Rabinovich et al. 2006*), that are well known in population biology and can explain the competition among three species: 

<img src="https://raw.githubusercontent.com/sazio/Transients/master/img/ML3D.png" width="200">



The variables x, y, z are the number of individuals for each population at time *t*, and <img src="https://render.githubusercontent.com/render/math?math=\rho_{ij}"> are competition coefficients measuring how much the *jth* species affects the growth rate of the *ith* species. The non-symmetry of the <img src="https://render.githubusercontent.com/render/math?math=\rho_{ij}"> guarantees the WLC behavior of the above mentioned dynamical system. The phase portrait of such a behavior is a heteroclinic contour (see below)


<img src="https://raw.githubusercontent.com/sazio/Transients/master/img/phase_portrait.png" width="400">



