[&larr; Back to TMD3.0 Main Page](../README.md)

# Getting Started with TMD
This page contains the code from the introductory video tutorial found on youtube. 

## 1. Installing TMD 
The first step toward getting started with TMD is to install the software on your computer. The installation process is described [here](doc/installing_tmd.md).

## 2. Getting tide model data 
After you've installed the TMD functions, you'll need some tide model data. [Click here for a list of all available tide model data](doc/tide_model_data.md). 

## 3. Getting help 
#### Online
The most up-to-date documentation can always be found on the main page of the [TMD GitHub repo](https://github.com/chadagreene/Tide-Model-Driver). There you can also post issues if you find any bugs in the code.  

#### In MATLAB

If you're in MATLAB and you want help with a specific function, you can get plain-text help in the Command Window by typing `help` followed by the name of the function. For example: 

```matlab
>> help tmd_predict
```
To access formatted documentation with lots of examples for any of the primary TMD functions, type `tmd` followed by the function name. For example:

```matlab
>> tmd tmd_predict
```
If you're not sure what function name you're looking for, just type `tmd` into the Command Window, and it will bring up a complete function list:

```matlab
>> tmd 
```

## 4. Predicting tides 
#### Single-location time series 

#### Drift-track time series 

#### Map-view snapshot 

## 5. Accessing tide model data 
#### Using `tmd_interp` to interpolate water column thickness

#### Using `tmd_data` to access gridded data


# Author Info
This page was written by [Chad A. Greene](https://www.chadagreene.com), January 2023.