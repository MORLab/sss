sss
====

A MATLAB toolbox for large-scale dynamical systems in state-space described by sparse matrices. It extends the capability of the Control System Toolbox by defining new Dynamic System Objects (sss) and analysis methods that exploit sparsity.

Working with large-scale dynamical systems has never been easier. All you need to do is add an "s" in the system definition

``sys = sss(A,B,C,D,E)``

and use all your favorite functions such as
`` bode(sys), step(sys), norm(sys), sys = sys1 - sys2, ...``
and many more!

For more information, type `doc` in the command window or visit http://www.rt.mw.tum.de/?sss. Check out also our demo by typing `sss_gettingStarted` in the command window.

***
*Programmed with:* MATLAB R2015b

*Tested on:* MATLAB R2014b, R2015b, R2016b (both Windows 7 and Ubuntu 16.04.1 LTS)

*Some functions require:* Control System Toolbox

> Note: Sign up for our newsletter under https://lists.lrz.de/mailman/listinfo/morlab to stay up to date.

***
Copyright
----------
This toolbox is developed by [MORLab](https://www.rt.mw.tum.de/?morlab), the model reduction lab at the [Chair of Automatic Control](https://www.rt.mw.tum.de/en/home/) in collaboration with the [Chair of Thermofluid Dynamics](http://www.tfd.mw.tum.de/index.php?id=5&L=1).

***
Acknowledgements
-----------------
The developing team is thankful to all the research assistants and students at [MORLab](https://www.rt.mw.tum.de/?morlab) that have contributed at creating and developing the sss class since 2008.

The team of [Morembs](http://www.itm.uni-stuttgart.de/research/model_reduction/MOREMBS_en.php), a model reduction software for elastic multibody systems, is sincerely acknowledged for the support in the automated generation of the documentation for the toolbox.

***
Developing guidelines
----------------------

We hope that you enjoy the toolbox and would like to contribute by extending its capability.
To make sure that the developing does not get out of hand, we prepared a few guidelines that we ask you to follow.


### Folder structure
The folder structure of the toolbox is as follows
- **sss** (main folder)
	- **src** (source code)
		- **@sss** (sss class and methods)
		- **+sssFunc** (functions used by different classes)
		- **extras** (additional functions)
		- **sim** (simulation)
    - **benchmarks**
    - **doc**
	- **demos**
	- **test**

### Documentation
To automatically generate the documentation for the toolbox from the function headers, type `publishDoc('sss')` in the command window. Make sure to format the function headers according to the ``headerTemplate.m`` provided. To publish the documentation for a single function, use syntax `publishFunction('function name')`.

### When should a method be added within the ``sss`` constructor or within the @sss folder?

Following functions/methods should be included in the class definition (constructor)
- set and get
- everthing that is small enough to fit in 5 lines (e.g. one if else end section plus assigments)
