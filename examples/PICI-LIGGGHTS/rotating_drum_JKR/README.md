# Using the JKR Model in PICI-LIGGGHTS

The JKR Contact Model is Cohesive Contact Model with the intention of replacing of the 'sjkr' (or simplified JKR) model in the standard LIGGGHTS library.

The implementation of this model was developed by Eidevåg, T., Abrahamsson, P., Eng, M., & Rasmuson, A. [1] based on the 2009 article "Discrete-element modeling of particulate aerosol flow" by J.S. Marshall [2].

## Liggghts simulation setup:

Simulation set up requires 2 main changes from a typical 'Hertz' set up:

-1. Use the 'gran model jkr' pair style
-2. Use two alternative 'fix property/global(s)'

-----------------------------------------------------------
Pair Style

1.     
Use your normal simulation file and change the 'pair_style' to include the JKR model:

'pair_style     gran model jkr'

-----------------------------------------------------------
Property/global(s)

2. 
2.1 Setting the Work of Adhesion:

    'fix         {fix_variable} all property/global workOfAdhesion peratomtypepair {mat_property}'

The work of adhesion is synomous with the 'Surface Energy', with units of J/m^2.

!! This property essentially defines how 'sticky' the particles are.
 
2.2 Setting the resolutionJKR scalar:

    'fix        {fix_variable} all property/global resolutionJKR scalar {mat_property}' 

!! The resloutionJKR scalar variable defines the "resolution" for a look-up table that is calculated at the start of a simulation. The resultant, speeds up the simulation by interpolating the contact force. 

The typical value for simulations is 1e-4, and Tobias suggests it should be sufficient for mostsimulations [3].

If the user is presented with: "Simulation outside LUT table, increase the resolution for speed! delta_delta_c = " see below:

Lowering the resolutionJKR gives a better resolution but lower range of displacement values.

What the error means: 
The value is outside the given look-up table will then solve the fourth-order polynomial for each iteration, significantly slowing the simulation. This typically happens in 'ill-conditioned' cases, therefore, the user is suggested to use a lower value to increase simulation speed.

For more information on this issue please see the original blog post: See Source [3]  

-----------------------------------------------------------
Key:

{fix_variables}   = the user designated string variable
{mat_property}    = the desired material property

-----------------------------------------------------------

## Credits
https://github.com/eidevag/LIGGGHTS-PUBLIC-JKR

## Sources:
[1] Eidevåg T, Abrahamsson P, Eng M, Rasmuson A. Modeling of dry snow adhesion during normal impact with surfaces. Powder Technology. 2020 Feb 1;361:1081-92.

[2] Marshall JS. Discrete-element modeling of particulate aerosol flows. Journal of Computational Physics. 2009 Mar 20;228(5):1541-61.

[3] Tobias Eidevåg, CFDEM Forum, LIGGGHTS PUBLIC V3.8.0 WITH JKR CONTACT MODEL (link) https://www.cfdem.com/forums/liggghts-public-v380-jkr-contact-model
