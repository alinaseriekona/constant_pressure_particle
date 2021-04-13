# Methane Pyrolysis and Particle Formation in an Ideal Constant Pressure Reactor
This sample probelm is inspired by reactor1.py example from the Cantera python library.<br>
Some basic particle formation steps are added to the code.<br>
By typing python3 reactor1.py --plot, the plots of particle volume fraction, CH4, A4R5C, and A4R5 mole fractions are plotted. <br>
A4R5C is the consumed mole fraction of A4R5 after particle formaion in each iteration. <br>
I am looking for a way to overwrite A4R5C on A4R5 in the gas object before going to next iteration. Is this possible?
