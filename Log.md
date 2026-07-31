# Log for my studies

## Sprint Week

- In order to get a random set of data in a range, use the np.random.uniform, which has the values for min and max, instead of random.random, which doesn't

- In order to do np arrays like js maps, best bet is to try a lambda function, and vectorize it beforehand if you need to. You don't always need to. I do wish they had a map function.

- Floating point number issues are usually why you see larger clifs in the graph, not a data problem. 

- Make sure to import all of the sub stuff you need for species package

- Brown dwarfs (especally in the L/T transition) don't always show cooler = redder. Around T I discovered that methane can form in large enough quantities to block some of the H filter I was using from the tutorial. 
  
- Reddening vector shows the direction that should occur due to to adding dust to the measured values.

- If you get a dwarf that isn't matching the reddening factor you're tracking, it's probably not that factor that's causing the issue. 

- Best place to check the filter names is https://svo2.cab.inta-csic.es/theory/fps/index.php Don't expect them to just be their name.

- There's a difference between the filter value of the filter only, and the whole system. The whole system contains information about QE, throughput, as well as the filter values. 

-  Filter values come back as a 2d array of wavelengths and transmissions

- Working through repro-ing, it's time to start learning about atmo, and how it impacts brown dwarfs and the measurements. 

- Adding new variables in models increases the time to calclate the various grid points, so adding a new variable instead of using a constant requires a significant scientific reason to do so. 

- Another key point is that there isn't a magic function that can be called to calculate the perfect values in a downloadable model. To get around this scientists check certain points, and then have functions interpolate around them. Is there a function that could make that work? Yes. But running that function is for super-computers, not for python. 

- K_zz stands for vertical atmospheric mixing. A high mixing value means the atmo has a larger varying chem, and can produce unexpected spectrum. 
  
- Download a full model causes a huge memory usage. Not really what you want to do, so you need to set a wavel_range and resolution, so that it cuts the amount of data that you need to work with. For what I'm doing, I don't need the whole thing. 

- I've noticed from clauding around that the elf-owl models don't go as high as the JWST filter. I need to go through and figure out how the paper dealt with that kind of infomration. 

- A good tutorial for finding synthetic photometry lives here: https://species.readthedocs.io/en/latest/tutorials/model_spectra.html It's not typesafe, but it works

- To get flux, you want to use the ReadModel.get_flux functions to pull that out. Specify which filter value so you can do a full color - color fit. 

- To get Mag, you want to use the same class but get_magnitude

- Sigh... time to remember that abs mag and apparent mag are what would be if the object was viewed from 10 pc away, which turns out to be the same if you specify 10 pc. Huzzah, not crazy for seeing the same number twice. 

- We're running into an issue where the model specified in the paper don't go through the final F1550C filter's effective wavelength. elf-owl ends at 15um, while the F1550C's effective range is (np.float64(14.912), np.float64(16.219)), so barely hitting that lower end, but not getting anywhere close the peak of that filter. I get an error on species when attempt to load the model for those wavelengths. 

- This might be a good example where they're doing some sort of custom model entry, or I'm missing a species param to force it to extrapolate, but I'm not seeing it on my end. I've also checked the original elf-owl paper (https://arxiv.org/pdf/2402.00756) to verify that the wavelength values reported by species are correct. From what I can tell, they are. 

- But looking at the paper Michelson et al 2025, it says they used the library to get that info. Looks like a good time to reach out to Zac / the paper authors to see how they got those values. 


- In order to unblock myself, I'm going to try and create the graph with simulated photometry using the sonora-diamondback-highres model. It's smaller, and doesn't have the same, or as many params, so you won't be able to get the same resolution, but at least we'll have enough for the color/color diagram


- After doing some looking, I was able to determine that the diamondback model isn't a good fit because it doesn't focus on the cold end of the spectrum. It's not really what I'm looking for here, because it doesn't fit the science I'm testing. 


## 2026-07-27

### Staring the week
- We've got the graph as far as I can go without additional information. 
- I've reached out to Ashley Messier with a request to get the galaxy information and how to understand how they got the final JWST filter
- It's time to start jumping on some of the technical learnings that I'll need as well as do some galaxy exploration. Maybe I can find the data in other places. 

### What's blocked, and on whom
- Blocked on getting the filter data, as well as the Galaxy data. 

### What I learned 
- One of the issues that I ran into while working on the graph is that we don't have data for some of the points. So we have a grid, but not every point on the grid is filled. 
- To combat that, we look at the data and only sample from the points we know we have. It makes the graph look segmented, but we can play around with it after I get some of the other data. 

### Next
- Start on the stats + astro work that I generated to help get me up to speed. 

## 2026-07-29

Heard back from Ashley! She was very helpful and gave some advice on how to get the galaxy models, as well as some pointers for where my exploration can go later in the fall. 

Basically, she's going to reach out to the PI to see how they were able to get the F1550C's filter data to work with the elf-owl models, but also gave me some excellent ideas on how to extend my work by verifying the spread factor on a lot of the new models. Facinating stuff, but stuff that I'll need to check out after I do some replication. 

The specific model she's referening is the brand new exo-remk26 model that just came out this year. The data can be found here: https://lesia.obspm.fr/exorem/YGP_grids/Exo-REMk26/High_Res_grid_2026/R200k_cloudyfsed_2026/

I've also been working on making the galaxy graph, and I've realized that it's often helpful to look at the original author's website instead of a mirror, as the mirror doesn't have to have the complete data set. I ran into this with the SWIRE subset, and the original author's (Polletta) distro. 

I did some research into how Vega can affect the magnitude readings that come back from species. Thankfully this filters out when doing a color-color analysis, as everything is being affected by the same Vega drops in MIR range. However, when you redshift too far, you might run into issues with magnitude adjustments. 

Also dug into the difference between F_λ and νF_ν, and why astronomers like the log log adjusted graphs in νF_ν, as it shows it adjusted for the log scale. 

## 2026-07-30
I wanted to see if we could attempt to extrapolate the elf-owl models in the colder temps (y models), to get the same (or approx the same data that Michleson was showing). In order to check this, we're going to check by plotting a power law to data we do have, and see how far we get in the residuals. 

Based on my research, we can do the following: logF=logA−αlogλ.

One issue that I'm running into is the time it takes to load the spectra of the 200 planets that I auto generate. These will be critical as we move towards doing larger explorations. I'm going to have claude dig into the internals of the species library to see if we can shorten the load times. 

Success on the claude front. Turns out you can load the data directly from the .hdf5 database, without going through the model.getData function. Implementation lives here: helpers/generate_planet_list.py -> get_planet_spectra. 

In order to make sure the implementation is correct, I also made a second version of the call in the helpers/generate_planet_list.py -> get_planet_spectra_via_species. The outputs of those checks will show on the notebook. 

The results showed no residual between the getData implementation and the db implementation on the t models, however, there is a slight difference in the y models, with a max abs diff of 1.962e-34, with a relative dif of 3.813e-16. Deviations in these degrees don't consititute an issue for the work we're doing here, as it's essentally a rounding issue in floating point.

Claude helped with those particular functions, as I was unaware of the internals of the .hdf5 file. However, it did speed up the retrieval time for those functions. 

To check the residuals, we're going to take the wavelengths and fluxes in the 12-14 um range, and calcuate a power law, and attempt to extrapolate through the 14 - 14.9um data that we already have. This will allow us to check to see if it's a viable strategy moving forward. 

In generating the data with the elf-owl models, we got the following data back. 

| model | p16 residual (%) | median residual (%) | p84 residual (%) |
| -----| ------| ----- | ----|
| sonora-elfowl-t | -14.0 | -0.7 | 22.5 |
| sonora-elfowl-y | -6.7 | 11.7 | 49.8 |

Here's what that tells us. 
1. For the t model, the median residual is about even, meaning that it's symetric in it's reporting, and could be tentativly used to further extrapolate, but only on the base line (no spectra features)
2. The y model shows a different story. It's median is >10x the size of the t model, and it gets worse the colder the dwarf gets. With a median of 11.7%, its fair to justify that the power law is attempting to extrapolate more flux than is actually being predicted. That means even for a baseline, extrapolation isn't viable.

To dig a bit deeper, we can check out the sonora-bobcat model, which has a temp range of 200-2400K and a 0.61-17um wavelength range, which covers the F1550C filter completely. We can do the same extrapolation, but in the filter's range to verify our assumption.  

We can re-use our functions for the previous generation and breakdown, which should speed up our ttd (time to data (trademark pending))

The code to do that extrapolation is located in michelson-bobcat-extrapolation-check.ipynb

The results validate our assumption. Here's the table binned to the specific ranges. This is taking a power law configured at the 12um to 14um range, and then doing the extrapolation / residual check on the bobcat model.

| wavelength range | median residual % |
|--|--|
| 14.0 - 14.9um | 6.1% |
| 14.9 - 15.5um | 5.2% |
| 15.5 - 16.0um | 25.2% |
| 16.0 - 16.6um | 14.6% |
| 14.912 - 16.219 (F1550C's coverage) | 10.7% |

This validates our use of the bobcat model as the 14-14.9um median is 6.1% which is bracketed by the elfowl t and y models.

All of these residuals are positive leaning, which validates our conclusion that extrapolation would cause an over approx of flux, while not accounting for spectral features that may help break apart the brown dwarf population with the red-shifted galaxy population. 

## 2026-07-31

Went back to the elf-owl residual data and binned it by Teff instead of wavelength, showcasing that the colder temps are causing issues with the extrapolation (which is the whole point of the paper)

**sonora-elfowl-t**

| teff (K) | n_planets | p16 residual (%) | median residual (%) | p84 residual (%) | p84 \|residual\| (%) |
|---|---|---|---|---|---|
| 575.0 | 21 | -14.9 | 1.3 | 36.8 | 36.8 |
| 600.0 | 13 | -13.0 | 0.3 | 26.6 | 31.5 |
| 650.0 | 16 | -13.5 | -0.9 | 24.8 | 28.6 |
| 700.0 | 13 | -14.2 | -0.5 | 25.9 | 28.0 |
| 750.0 | 17 | -19.5 | -3.1 | 19.3 | 29.7 |
| 800.0 | 16 | -13.8 | 0.4 | 24.6 | 26.3 |
| 850.0 | 20 | -14.3 | -0.6 | 22.4 | 26.6 |
| 900.0 | 22 | -18.5 | -4.3 | 17.7 | 26.0 |
| 950.0 | 9 | -17.7 | -3.5 | 16.2 | 23.5 |
| 1000.0 | 19 | -14.8 | -2.4 | 17.0 | 22.9 |
| 1100.0 | 13 | -14.4 | -3.1 | 14.3 | 21.2 |
| 1200.0 | 21 | -14.7 | -3.2 | 15.3 | 20.9 |

**sonora-elfowl-y**

| teff (K) | n_planets | p16 residual (%) | median residual (%) | p84 residual (%) | p84 \|residual\| (%) |
|---|---|---|---|---|---|
| 275.0 | 9 | 17.2 | 41.4 | 83.7 | 83.7 |
| 300.0 | 15 | 7.3 | 25.7 | 64.7 | 64.7 |
| 325.0 | 24 | 2.6 | 21.2 | 56.1 | 56.1 |
| 350.0 | 19 | 1.1 | 18.0 | 45.4 | 45.4 |
| 375.0 | 12 | -2.7 | 11.6 | 40.7 | 40.7 |
| 400.0 | 17 | -11.7 | 9.2 | 44.1 | 44.1 |
| 425.0 | 14 | -5.9 | 10.4 | 34.3 | 35.4 |
| 450.0 | 15 | -11.5 | 6.2 | 37.9 | 37.9 |
| 475.0 | 19 | -13.2 | 3.8 | 34.5 | 36.2 |
| 500.0 | 17 | -11.4 | 5.2 | 31.6 | 32.5 |
| 525.0 | 20 | -13.9 | 0.9 | 32.8 | 34.0 |
| 550.0 | 19 | -14.5 | 2.2 | 36.7 | 36.7 |

