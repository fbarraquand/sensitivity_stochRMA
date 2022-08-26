# Investigating structural sensitivity in the stochastic Rosenzweig-MacArthur with environmental noise

Companion code to *No sensitivity to functional forms in the strongly forced, continuous-time stochastic Rosenzweig-MacArthur model*, Frédéric Barraquand. 

* ``main_figures`` reproduce Figs. 1, 2 and 3. 
* ``largeK`` (up to $K=15$), ``sigma_E=0.1`` (smaller noise) and ``noisyFR`` (add noise on the functional response) those of Appendices 1 and 2. 
* ``epsilon=0.3`` further reduces the conversion efficiency $\epsilon$ (not done in Fussman and Blasius 2005). This is meant to examine structural sensitivity for a slightly more realistic model. The results are very similar to $\epsilon = 1$. This code has been rewritten in C++ in ``SDE_CPP``. 



