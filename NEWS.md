# AHPtools 0.2.0

Changes:

* New function revisedConsistency() added. Takes a PCM as the parameter

* `createLogicalPCM()` which was in a separate package by the same name has
  now been included in this version

* All functions in the package that take a PCM as the input parameter
  have been modified to alternatively take a vector of the upper triangular
  elements of the PCM. A second parameter to these functions, 
  typePCM (TRUE/FALSE), is used to indicate if the input is a PCM or a vector
  of the corresponding upper triangular elements.
