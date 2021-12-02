# RR-Lyrae-Research
Curve fitting analysis for candidates of a RR Lyrae - binary star system

This was my final addition to an ongoing research project at Grove City College to test the binarity of RR Lyrae stars.

## Theoretical Model
RR Lyrae are variable stars that having pulsating brightness levels on the order of a few hours to a few days. Binary systems consist of two stars that orbit around each other with a shared center of mass and whose orbits can have a period from mere minutes to hundreds of years. In the case of this research, we primarily were focused on binary systems with periods estimated on the order of 10-100 years. There have been confirmed cases of binary systems in which one of the stars is an RR Lyrae. These cases have been observed using the Light Travel Time Effect (LiTE) in which the times when the RR Lyrae reaches peak brightness, also referred to as times of maxima, are recorded over a period of a few years. If the RR Lyrae is a member of a binary system, then when it is in the apogee of its orbit relative to Earth, the time it takes the light to reach us will be longer than when it is at its perigee. We can plot these times of maxima and see a sinusoidal curve as the star completes one orbit. If we can fit an equation to this curve, we can extract the period of this binary system and extrapolate this to other characteristics of the stars such as mass and luminosity.

## Mathematical Approach
Research has been done in this field prior and is referenced in the 'papers' section of this repository. Specifically, the paper written by *J.-L. Li, et. all*, describes this curve of light maxima with the equation below. 

O - C = &Delta;T<sub>0</sub> + &Delta;P<sub>0</sub> * E + <sup>&beta;</sup>&frasl;<sub>2</sub> E<sup>2</sup> + &tau;

The goal of our research was to primarily find a method to fit the curve from data to this equation, henceforth referred to as the LiTE equation. Being as this is the summary for the code and not a complete research paper, users wishing to gain a greater understanding of the physics behind this phenomenon may reference the aforementioned papers for additional information.

## Curve Fitting
The LiTE equation (including further extrapolations) contains seven variables that must be solved to generate a perfect fit. When solving for particular variables, we compared our calculated results with those calculated in *J.-L. Li, et. all*. These seven variables have been broken into three categories for simplicity of explanation. They are as follows:
1. The linear portion
2. The parabolic portion
3. The binary portion

The **linear portion** contains the variables &Delta;T<sub>0</sub> = and &Delta;P<sub>0</sub> which are the initial epoch of the data sample (can arbitrarily be set to zero) and the pulsation period of the variable star. These two dominate the data initially and cause the O-C diagram (observed data minus calculated) to take on a form similar to
O - C = &Delta;P<sub>0</sub>E + &Delta;T<sub>0</sub>. 

![image](https://user-images.githubusercontent.com/38231105/144331624-67a7af48-1e86-4360-8fbb-e92014aac2d3.png)

As you may guess, a simple curve fitting algorithm method can be used to find these two variables. The find these values, we implemented a least-squares method to solve. More information on this type of curve-fitting method is widely available. <a href=https://en.wikipedia.org/wiki/Least_squares>(wiki link)</a>

Once the initial two values are calculated, we subtracted them from the data, as the term observed-calculated infers, and are presented with a filtered data results that looks like this.

![image](https://user-images.githubusercontent.com/38231105/144332174-406f7b5f-776f-4acf-a818-f068abd42149.png)

This is what we refer to as the **parabolic portion**. Though it is not as easily visible as the linearity of the initial data was, this data does have a slight positive curvature overall which can be modeled in the form &beta;E<sup>2</sup> + A E + C where &beta; is the same as in the LiTE equation, E is the epoch, and the rest are extraneous variables required to generate a basic parabola. This portion too was easily solved using the least-squares method and &beta; was solved for.

The **binary portion** contains the remaining four variables and was the primary hurdle of this research. Leaving off the linear and parabolic portions of the data, the O - C diagram looks like this.

![image](https://user-images.githubusercontent.com/38231105/144333497-20faa85d-0606-4f0c-99cc-f7fbedf06642.png)

As indicated by the blue line, this data takes a sinusoidal form that could be entirely a result of the RR Lyrae being one member of a binary system. The equation for this portion of the graph, however, is not a simple sinusoid, but rather is defined as follows.

&tau; = A\[&radic;(1-e<sup>2</sup>) sinE<sup>\*</sup> cos&omega; + cosE<sup>\*</sup> sin&omega;\]

Because this equation has four variables for which we need to solve, the solution space is much larger than before. We initially attempted a least-squares method as with the other two portions, but this did not yield the results that were consistent with *J.-L. Li, et. all.* Instead, we began to use a <a href=https://en.wikipedia.org/wiki/Radial_basis_function>radial basis function optimization</a> algorithm (RBFopt) to approximate our results. This is the function that can be seen throughout the code. In short, 
RBFopt shrinks the possible solution space but performs more complex steps to approximate the solutions, so even though it may be more accurate than the least-squares method, it takes much longer to run.

By the end of the semester, we still had not created a successful implementation of RBFopt and I passed down the project to an underclassman who is carrying on the project now. The next approaches I suggested were to explore other mathematical models that could represent the LiTE curves. So, while this code is ultimately incomplete, it does demonstrate a systematic approach to breaking down an equation with a large solution space and using different methods to appoximate solutions to a curve.


## Additional notes on code and data
The data we used originated from the <a href=http://rr-lyr.irap.omp.eu/dbrr/rrdb-V2.0_08.3.php?XX+And&>GEOS RR Lyrae database</a>. A porttion of my time was spent using the <a href=https://pandas.pydata.org>pandas library</a> and Microsoft Excel to clean the data from this database. I also utilized <a href=https://matplotlib.org>matplotlib</a> to graph our data and the curves that we generated.
