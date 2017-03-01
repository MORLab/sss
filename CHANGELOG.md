sss - Changelog
================

A list of (major) changes between releases. Sometimes we add also changes to come to our **roadmap**.
***

Roadmap (changes to come)
-------------------------
- pss: parametric models
- second order models
***

v1.03 [02 February 2017]
-----------------------
### Changes
- LYAPCHOL
	* Changed nomenclature to be more consistent with the literature.
	* Changed definition of (low-rank) Cholesky factors from upper triangular to lower triangular to be better suited for the low-rank case.
	* Option `type` was renamed to `method`. Optional values now include `'auto', 'adi' and 'hammarling'`.
	* Added options `forceOrder` and `maxiter`.

-	RESIDUE
	* Computing only the first 6 residues corresponding to the eigenvalues with smallest magnitude, instead of all residues and in previous versions, which is too expensive in large-scale setting.
	* Added option `eigs` to specify the type of eigenvalues to be computed `'sm','lm',...`
	* Added option `nEigs` to specity the number of pole/residues to be computed.

- DECAYTIME
	* Adapted to used the new version of `residue`

- SSS
	* Added property `isSym` which checks for symmetry of A and E
	* properties such as `isSiso, isMimo, isBig, isDae, isDescriptor` are no more hidden.

- SECOND2FIRST
	* Added system matrices A,B,C,D,E as possible outputs

### Bugfixes
	- MINUS
		* changed the function to be able to compute the difference of sss and ss objects as well.
	- ZEROS, ZPK
		* usage of sys.p and sys.m corrected
	- LOADSSS
		* defining zero matrices with `spalloc` instead of `zeros` to reduce memory consumption.

### Third-party
- MMESS 1.0.1
	* The integration of the MMESS toolbox was changed to be able to use MESS without modifications from our side. This should allow you to update MMESS independently of the sssMOR toolbox. However, compatibility can be guaranteed only with the version we used for development (in this release: MMESS 1.0.1). In addition, solveLse needs to be added to the usfs.
***

v1.02 [05 October 2016]
-----------------------

### Changes
- NORM
	* Adding a "stabcheck" option to be able to avoid the (sometimes expensive) stability check before the norm computation
***

v1.01 [16 September 2016]
-------------------------
### Changes
- SOLVELSE
	* **new function: solve linear system of equations**

- LYAPCHOL
	* **new function: solve Lyapunov equation**

- PZMAP
	* plot only the most important poles and zeros

- ZPK
	* **new function: create zpk objekt of the most important poles and zeros**

- POLES
	* **new function: compute the most important poles**

- ZEROS
	* **new function: compute the most important zeros**

- DCGAIN
	* **new function: compute dc gain of lti system**
***

v1.00 - first Release [16 November 2015]
-----------------------------------------
