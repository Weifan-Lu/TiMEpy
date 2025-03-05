## Adjusting the Sampling Rate in TidalStrain.2

The TidalStrain.2 project allows you to customize the sampling rate for tidal strain calculations by modifying parameters across several key files. This section guides you through the process step-by-step, ensuring consistency and accuracy in your adjustments.

### Before Changing the Sampling Rate, Pay Attention to:
**TidalStrain.2.sh**: The main script file where the `nevent` parameter needs to be modified.  
**event.in**: The input file containing event data, with the number of rows matching the `nevent` value.


## Overview of Relevant Files

To change the sampling rate, you’ll need to update the following files:

-  **earthtide.input5.base**: Controls the solid Earth tide calculations.
-  **ocean_loading.in.base**: Manages the ocean tide loading effects.
-  **sum_both_list.f**: A Fortran script that combines results from the above calculations.

## Step-by-Step Guide to Changing the Sampling Rate

Follow these steps to adjust the sampling rate to your desired interval:

1. **Modify | `earthtide.input5.base`**

```base
    JA,JB,JC,JD,JE :initial Year,Month,Day,Hour,Minute
    961,0.05 :number of calc., calc. interval(hours)
    PH,RM,0.0,DP :lat.,lon.,height(m),depth(km)
    56.0,0.D0 :ET-UT(sec), diff. from UT(hour)
    6 :Kind(1:grv,2:NStlt,3:EWtlt,4:vol,5:sea,6:strain-tensor)
```
NN: Number of calculations. For example, if you would like to set the sampling interval at 3 minutes for 2 days, you should set 961 (= 60*24*2 / 3 +1). If 6 minutes for 1 day, 241 (= 60*24 / 6 +1). If 15 minutes for 3660 days (about 10 years), 351361 (= 60*24*3660 / 15 +1).
ST: Step or interval of calculations in hours. Note that we recommend 0.05 (= 3 minutes) as minimum because the unit of the sampling interval of ocean tide loading effects combined later is integer minutes.

2. **Modify | `ocean_loading.in.base`**

```base
    *********************[ Mandatory Cards ]**********************
    STAPOSD eventID , longitude , latitude , depth , 0.0
    WAVE ALL
    KIND ST
    **********************[ Option Cards ]************************
    PREDICT 1,starttime,endtime,3
    PREFMT 5,1
    GREENF 1
    MESH3 ON
    MESH4 ON
    FULLMESH ON
    UNIT6 ocean_loading.log
    UNIT20 ocean_loading.out
    END
```

If you would like to set the sampling interval to 15 minutes, you should change 3 in “PREDICT 1,starttime,endtime,3” to 15.


3. **Modify | `sum_both_list.f`**

If you would like to change the sampling interval and duration, you can change red values of “tinterval=3.0” and “do i=1,961” in sum_both_list.f.

For example, every 6 minutes for 1 day: 
```
tinterval=6.0, do i=1,241.
```

1 minutes for 28 days: 
```
tinterval=1.0, do i=1,40321.
```

Of course, you must change the following red values properly.

- Solid tide: NN and ST (961,0.05) in “earthtide.input5.base (see 1. **Modify | `earthtide.input5.base`**)”
- Ocean tide: PREDICT 1,starttime,endtime,3 in “ocean_loading.in.base (see 2. **Modify | `ocean_loading.in.base`**)”
