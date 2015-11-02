# sss
A MATLAB toolbox for large-scale dynamical systems in state-space described by sparse matrices. It extends the capability of the Control System Toolbox by defining new Dynamic System Objects (sss) and analysis methods that exploit sparsity.

Working with large-scale dynamical systems has never been easier. All you need to do is add an "s" in the system definition

``sys = sss(A,B,C,D,E)``

and use all your favorite functions such as

``bode(sys),
step(sys),
norm(sys),
sys = sys1 - sys2,
...``

and many more! 

Check out our documentation (type ``doc`` in the command window) and demos to get started.

**COPYRIGHT**
This toolbox is developed by [MORLab](http://www.rt.mw.tum.de/en/research/fields-of-research/model-order-reduction/), the model reduction lab at the [Chair of Automatic Control](www.rt.mw.tum.de/en) in collaboration with the [Chair of Thermofluid Dynamics](http://www.tfd.mw.tum.de/index.php?id=5&L=1).

**TERMS AND CONDITIONS**

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either GPLv2 or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

**Programmed and tested with: MATLAB R2015b**

*Some functions require the Dynamic Systems Toolbox*

**ACKNOWLEDGEMENTS**
The developing team is thankful to all the research assistants and students at [MORLab](http://www.rt.mw.tum.de/en/research/fields-of-research/model-order-reduction/) that have contributed at creating and developing the sss class since 2008. 

The team of [Morembs](http://www.itm.uni-stuttgart.de/research/model_reduction/MOREMBS_en.php), a model reduction software for elastic multibody systems, is sincerely acknowledged for the support in the automated generation of the documentation for the toolbox.

# Developing guidelines
We hope that you enjoy the the toolbox and would like to contribute by extending its capability. 
To make sure that the developing does not get out of hand, we prepared a few guidelines that we ask you to follow. 


### Folder structure
The folder structure of the toolbox is as follows
- **sss** (main folder)
	- **src** (source code)
		- **@sss** (sss class and methods)
		- **extras** (additional functions)
	- **demos**
	- **test**

### Documentation
The documentation for the toolbox is automatically generated from the function headers. Make sure to format the header according to the ``headerTemplate.m`` provided.

The functions to publish the documentation are not included in the git repository as they belong to a third party.

### When should a method be added within the ``sss`` constructor or within the @sss folder?

Following functions/methods shuld be included in the class definition (constructor)
- set and get 
- everthing that is small enough to fit in 5 lines (e.g. one if else end section plus assigments)

