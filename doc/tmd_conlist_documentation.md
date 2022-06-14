[&larr; Back to TMD3.0 Main Page](../README.md)

# `tmd_conlist` documentation

`tmd_conlist` returns a list of tidal constituents in a TMD3.0 compatible
consolidated NetCDF tide model file. 

## Syntax

 conList = tmd_conlist(filename)

## Description 

`conList = tmd_conlist(filename)` returns a cell array of tidal
constituents in the specified model file. 

## Example 
Get a list of constituents in the updated CATS2008 model: 

```matlab
>> conList = tmd_conlist('CATS2008_update_2022-06-11.nc')
conList =
  1×10 cell array
    {'m2'}    {'s2'}    {'n2'}    {'k2'}    {'k1'}    {'o1'}    {'p1'}    {'q1'}    {'mf'}    {'mm'}

```
## Author Info 
The `tmd_conlist` function and its documentation were written by Chad A.
Greene, June 2022. 