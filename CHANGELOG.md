sss - Changelog
================

A list of (major) changes between releases, as well as information about MATLAB versions used and toolbox dependencies. Sometimes we add also changes to come to our **roadmap**.
***

Roadmap (changes to come)
-------------------------
- second order models
***

v2.00 [06 September 2017]
-----------------------

|                 |                                                                                        |
|:----------------|:---------------------------------------------------------------------------------------|
| Dependencies    | Control System T., System Identification T. Optimization T., M-MESS, uniquetol, mmread |
| Programmed with | MATLAB R2015b, R2016b,                                                                 |
| Tested with     | MATLAB R2014a, R2015b, R2016b, R2017a                                                  |
| on              | Windows 7                                                                              |

### New Features
- DOCUMENTATION
  * added documentation files to a folder doc/ inside sss. In this way, everybody can profit from the doc documentation of sss, not just those who download the release version.
  * added p-Functions publishDoc.p and publishFunction.p to be able to update the documentation whenever headers are changed.
- SIMINIT, SIMUPDATE
  * added auxiliary functions in the folder sim/ that are used for the initialization and update of the simulation functions. Common computations are thus performed in these new functions.
- @SSS, +SSSFUNC
  * certain sss-functions (like e.g. `'decayTime', `'diag', `'eigs', `'lyapchol', etc.) can be now also passed with ssRed objects. This has the advantage that both full order models (sss objects) as well as reduced order models (ssRed objects) can be analyzed with our sss functions.

### Changes
- SSS
  * updated to v2.00
  * SSSMOR is no more a submodule of sss, hence it is not anymore in src/ folder. Instead, sss and sssMOR are now distributed independently and have to be both in the path for sssMOR to work.
- LOADSSS
  * loadSss is deprecated and will be removed in later releases of sss. Use |sss(fname)| instead.

v1.03 [02 February 2017]
-----------------------
|                 |                                                                      |
|:----------------|:---------------------------------------------------------------------|
| Dependencies    | Control System T., System Identification T., Optimization T., M-MESS |
| Programmed with | MATLAB R2015a, R2015b                                                |
| Tested with     | MATLAB R2014b, R2015b                                                |
| on              | Windows 7, Ubuntu 16.04.1 LTS                                        |
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
|                 |                                                              |
|:----------------|:-------------------------------------------------------------|
| Dependencies    | Control System T., System Identification T., Optimization T. |
| Programmed with | MATLAB R2015a, R2015b                                        |
| Tested with     | MATLAB R2014b, R2015b                                        |
| on              | Windows 7, Ubuntu 16.04.1 LTS                                |

### Changes
- NORM
	* Adding a "stabcheck" option to be able to avoid the (sometimes expensive) stability check before the norm computation
***

v1.01 [16 September 2016]
-------------------------
|                 |                                                              |
|:----------------|:-------------------------------------------------------------|
| Dependencies    | Control System T., System Identification T., Optimization T. |
| Programmed with | MATLAB R2015a, R2015b                                        |
| Tested with     | MATLAB R2014b, R2015b                                        |
| on              | Windows 7, Ubuntu 16.04.1 LTS                                |

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
|                 |                                                              |
|:----------------|:-------------------------------------------------------------|
| Dependencies    | Control System T., System Identification T., Optimization T. |
| Programmed with | MATLAB R2015a, R2015b                                        |
| Tested with     | MATLAB R2014b, R2015b                                        |
| on              | Windows 7, Ubuntu 16.04.1 LTS                                |
