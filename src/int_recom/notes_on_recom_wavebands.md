*** Preprocessor versus logical flags ***

The code for calculating spectral irradiance etc. has a number of switches. 
In the original MITgcm code version, most if not all of them were coded as
preprocessor statements. In the version for FESOM-REcoM only one of them
remains as a preprocessor statement, all others have become logical switches
are defined in recom_moduls.F90 and can be changed via namelist settings.

To me this makes complete sense. There are, however, two small caveats here:

- first, the compiler switch is sometimes called `RECOM_WAVEBANDS`, and 
  sometimes `__RECOM_WAVEBANDS`. The latter is the standard convention in FESOM
  and is also, how the switch is called in the file CMakeLists.txt  

- second, the switch to logical variables is not completely consistent.
  In recom_forcing.F90 we still find `#ifdef RECOM_CALC_REFLEC` although
  this should now be a logical variable, i.e. this should be a fortran-if

These two things should be corrected, and then the code should be re-compiled
and a testrun made.

*** Different options for using spectrally resolved light field in REcoM ***

The fortran logical switches defined within the RECOM_WAVEBANDS preprocessor
statement are (with their default values):

- RECOM_CDOM       = .true.: This switches on that an additional prognostic
  tracer is calculated. I am not sure whether this works. It should a) define
  the new tracer (so that e.g. an initial condition is read in), b) define parameters
  like a fraction of DOC production that goes into CDOM, and CDOM degradation, and 
  c) define terms in recom_sms.F90 that describe production and decay of CDOM. 
  Need to check! 

- RECOM_MARSHALL   = .false.: This switches on the mechanistic description of
  photoinhibition by Marshall and Geider. It also defines new tracers, namely the
  state of the D1 protein complex (degraded/intact) for each phytoplankton. I 
  don't know whether this works, but at the moment we will not use it. 

- RECOM_RADTRANS   = .false.: This switches between a simple attenuation
  calculation and the full radiative transfer code from Dutkiewicz. This we want
  to use!

- OASIM            = .false.: This switches beween reading in separate forcing
  files for the different wavebands, and the calculation of the waveband radiation
  from a 'typical' spectrum of PAR. 

- RECOM_CALC_ACDOM = .false.:

- RECOM_CALC_APART = .false.:

- RECOM_CALC_APHYT = .false.:

*** Some steps to remove problems with the code ***

When consistently using the preprocessor option __RECOM_WAVEBANDS parts of the code get
included that so far have not been tested, and it turns out that large parts of the
code do not properly work yet. Here are steps that needed to be done

recom_forcing.F90:

- replace `#ifdef RECOM_WAVEBANDS` with `#if defined (__RECOM_WAVEBANDS)`

- several variable definitions were encluosed in if statements. That does not work in fortran, like
  with preprocessor statements. In consequence these fields are now always defined, even if
  they are not needed

- several if statements were not working: e.g.: `if (RECOM_RADTRANS=.true.)` should be replaced with
  `if (RECOM_RADTRANS) then` ('=' is not a logical operator, and sometimes the 'then' was missing)

- in the call to recom_sms additional arguments are needed. That's fine. but in the call,
  some arguments have been enclosed in `if (enable_coccos) then`. That only works for
  a pre-processort statement, not within fortran!

- some MITgcm-like constants were used; `XXX .gt. 0 _d 0` needs to be replaced by `XXX .gt. 0.0d0`

- at one pooint there was an endif too many




