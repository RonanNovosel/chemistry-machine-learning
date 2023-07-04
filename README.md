# chemistry-machine-learning




This is my own personal GitHub page, which contains some of the machine learning, cheminformatics and computational chemistry code and data I use

This is purely for personal use, both as a student and a proffesional chemist

~ xoxoxo


The main projects currently are: 

- A simple test of a cheminformatics predictor for LogP using SMILES as the input. This is done mainly as a tutorial, but I will copy the code and expand the data set to reflect a wider spectrum of chemical space. Also, I would like to optimise the descriptors that are most relevant to the problem. I also would like to analyse and use a different machine learning technique, as well as computational chemistry to determine new, specific descriptors.

- A Chemical space viewer, allowing the view of different properties as related to descriptors such as MW, polar surface area and rotatable bond number, etc. The dots are then coloured on a 3D map to show the relevant properties. The value of this would be identifying unfilled (possible) areas of the chemical space to use in order to improve the accuracy of the model. This can be used in combination with automated computational chemistry calculations for harder-to-reach descriptor values.

- A machine learning predictor for the yield of Suzuki-Miyuara reactions. The data was exported from Reaxys and is used to analyse the SMARTS pattern of the reaction. It also takes into account temperature (degrees C), reaction time (hours), solvent, and catalyst. This will allow for easy prediction of the yield from a given theoretical reaction.
